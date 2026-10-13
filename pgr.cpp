// SPDX-License-Identifier: LGPL-2.1-or-later
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"
using namespace gmx;

class Pgr : public TrajectoryAnalysisModule {
    SelectionList sel_;
    std::string out_ = "pgr.xvg";
    double bin_=0.01, rmax_=2.0, bulkStart_=-1, bulkEnd_=-1;
    int samples_=200000, seed_=2026, frames_=0;
    std::vector<int> oxygen_;
    std::vector<uint64_t> counts_, hits_;
    std::vector<double> volumes_;
    uint64_t bulkCount_=0, bulkHits_=0, outside_=0;
    double bulkVolume_=0, bulkDensity_=0;
    std::mt19937_64 rng_;
    AnalysisNeighborhood nb_;
public:
    void initOptions(IOptionsContainer* o, TrajectoryAnalysisSettings* s) override {
        static const char* const help[]={
            "Proximal RDF: each selected water oxygen contributes once at its minimum",
            "periodic distance to ANY atom in the entire target selection.[PAR]",
            "Select exactly two groups with -sel: target first, water second.",
            "Water oxygen atoms are identified by atomic number 8, or by names starting",
            "with O when element data are absent (GRO). Select a water-only group.[PAR]",
            "Uniform random points in each frame's periodic box estimate proximal shell",
            "volumes. rho(r)=sum water counts / sum shell volumes; pg(r)=rho(r)/rho_bulk.",
            "The bulk interval defaults to [0.8*rmax,rmax). Set it explicitly after",
            "checking the density plateau. Unsampled bins are NaN, never zero.",
            "A companion OUTPUT_density.xvg contains density, summed volume, water counts",
            "and random hits. Increase -samples and check convergence with -seed.[PAR]",
            "Distances are atom-center distances in nm; atomic radii are not subtracted.",
            "Three-dimensional periodic boundary conditions are required."
        };
        s->setHelpText(help);
        s->setFlag(TrajectoryAnalysisSettings::efRequireTop);
        o->addOption(FileNameOption("o").filetype(OptionFileType::Plot).outputFile().store(&out_).defaultBasename("pgr").description("Normalized proximal RDF"));
        o->addOption(SelectionOption("sel").storeVector(&sel_).valueCount(2).required().onlyAtoms().onlyStatic().description("Target atoms first; water atoms second"));
        o->addOption(DoubleOption("bin").store(&bin_).description("Bin width (nm)"));
        o->addOption(DoubleOption("rmax").store(&rmax_).description("Maximum analyzed distance (nm), exclusive"));
        o->addOption(DoubleOption("bulk-start").store(&bulkStart_).description("Bulk interval lower bound (nm); -1 uses 0.8*rmax"));
        o->addOption(DoubleOption("bulk-end").store(&bulkEnd_).description("Bulk interval upper bound (nm); -1 uses rmax"));
        o->addOption(IntegerOption("samples").store(&samples_).description("Uniform volume samples per frame"));
        o->addOption(IntegerOption("seed").store(&seed_).description("Random seed for reproducibility"));
    }
    void initAnalysis(const TrajectoryAnalysisSettings&, const TopologyInformation& top) override {
        if(!std::isfinite(bin_) || !std::isfinite(rmax_) || bin_<=0 || rmax_<=0 || bin_>rmax_ || rmax_/bin_>1000000 || samples_<=0)
            GMX_THROW(InvalidInputError("Require 0 < bin <= rmax, at most 1,000,000 bins, and samples > 0."));
        if(bulkStart_==-1) bulkStart_=0.8*rmax_;
        if(bulkEnd_==-1) bulkEnd_=rmax_;
        if(!std::isfinite(bulkStart_) || !std::isfinite(bulkEnd_) || bulkStart_<0 || bulkStart_>=bulkEnd_ || bulkEnd_>rmax_)
            GMX_THROW(InvalidInputError("Require 0 <= bulk-start < bulk-end <= rmax."));
        const auto* atoms=top.atoms();
        std::set<int> target(sel_[0].atomIndices().begin(),sel_[0].atomIndices().end());
        std::set<int> residues;
        for(int i=0;i<sel_[1].posCount();++i) {
            int id=sel_[1].atomIndices()[i];
            const auto& a=atoms->atom[id];
            std::string name=*atoms->atomname[id];
            const bool isO=a.atomnumber==8 || (a.atomnumber==0 && !name.empty() && (name[0]=='O' || name[0]=='o'));
            if(!isO) continue;
            if(target.count(id)) GMX_THROW(InvalidInputError("Target and water oxygen selections overlap."));
            if(!residues.insert(a.resind).second) GMX_THROW(InvalidInputError("Multiple oxygen atoms in one selected residue: select a water-only group."));
            oxygen_.push_back(i);
        }
        if(oxygen_.empty()) GMX_THROW(InvalidInputError("No oxygen atoms found in the second selection."));
        size_t n=static_cast<size_t>(std::ceil(rmax_/bin_));
        counts_.assign(n,0); hits_.assign(n,0); volumes_.assign(n,0);
        // Expand cutoff slightly so an absent neighbor is never assigned to the last bin.
        nb_.setCutoff(static_cast<real>(rmax_+std::max(1e-5,1e-5*rmax_)));
        rng_.seed(static_cast<uint64_t>(seed_));
        fprintf(stderr,"pgr: %d target atoms, %zu water oxygens; bulk [%.6g, %.6g) nm; %d samples/frame\n",sel_[0].atomCount(),oxygen_.size(),bulkStart_,bulkEnd_,samples_);
    }
    void analyzeFrame(int, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData*) override {
        if(!pbc || pbc->pbcType!=PbcType::Xyz) GMX_THROW(InvalidInputError("pgr requires XYZ PBC; do not use -nopbc."));
        const auto& b=fr.box;
        double vol=std::abs(double(b[0][0])*(double(b[1][1])*b[2][2]-double(b[1][2])*b[2][1])-double(b[0][1])*(double(b[1][0])*b[2][2]-double(b[1][2])*b[2][0])+double(b[0][2])*(double(b[1][0])*b[2][1]-double(b[1][1])*b[2][0]));
        if(!std::isfinite(vol) || vol<=0) GMX_THROW(InvalidInputError("Invalid periodic box volume."));
        auto search=nb_.initSearch(pbc,sel_[0]);
        for(int i:oxygen_) {
            double d=search.minimumDistance(sel_[1].position(i).x());
            if(d<rmax_) ++counts_[std::min(size_t(d/bin_),counts_.size()-1)]; else ++outside_;
            if(d>=bulkStart_ && d<bulkEnd_) ++bulkCount_;
        }
        std::vector<uint64_t> h(hits_.size(),0);
        uint64_t bh=0;
        for(int i=0;i<samples_;++i) {
            double u[3];
            for(double& v:u) v=(rng_()>>11)*0x1.0p-53;
            rvec x;
            for(int k=0;k<3;++k) x[k]=u[0]*b[0][k]+u[1]*b[1][k]+u[2]*b[2][k];
            double d=search.minimumDistance(x);
            if(d<rmax_) ++h[std::min(size_t(d/bin_),h.size()-1)];
            if(d>=bulkStart_ && d<bulkEnd_) ++bh;
        }
        for(size_t k=0;k<h.size();++k) {hits_[k]+=h[k]; volumes_[k]+=vol*double(h[k])/samples_;}
        bulkHits_+=bh; bulkVolume_+=vol*double(bh)/samples_;
        ++frames_;
        if(frames_%10==0) fprintf(stderr,"pgr: processed %d frames\n",frames_);
    }
    void finishAnalysis(int) override {
        if(frames_==0 || bulkVolume_<=0 || bulkCount_==0)
            GMX_THROW(InvalidInputError("No usable bulk density: choose a populated bulk interval/rmax or increase samples. No normalized output written."));
        bulkDensity_=bulkCount_/bulkVolume_;
        fprintf(stderr,"pgr: %d frames, bulk density %.9g O/nm^3, bulk water count %llu, bulk MC hits %llu; %llu water observations beyond rmax\n",frames_,bulkDensity_,(unsigned long long)bulkCount_,(unsigned long long)bulkHits_,(unsigned long long)outside_);
        if(bulkHits_<10000) fprintf(stderr,"WARNING: fewer than 10000 bulk random hits; increase samples and check convergence.\n");
    }
    void writeOutput() override {
        std::string stats=out_;
        if(stats.size()>=4 && stats.substr(stats.size()-4)==".xvg") stats.resize(stats.size()-4);
        stats+="_density.xvg";
        std::ofstream o(out_), s(stats);
        if(!o || !s) GMX_THROW(FileIOError("Cannot open output files."));
        o<<std::setprecision(12); s<<std::setprecision(12);
        for(auto* f:{&o,&s}) *f<<"# pgr: minimum periodic atom-center distance, volume-corrected\n# frames "<<frames_<<" samples_per_frame "<<samples_<<" seed "<<seed_<<"\n# bulk_interval_nm "<<bulkStart_<<" "<<bulkEnd_<<" bulk_density_O_nm-3 "<<bulkDensity_<<"\n# bulk_water_count "<<bulkCount_<<" bulk_MC_hits "<<bulkHits_<<" bulk_volume_sum_nm3 "<<bulkVolume_<<"\n@ xaxis label \"Minimum distance (nm)\"\n";
        o<<"@ title \"Proximal radial distribution function\"\n@ yaxis label \"pg(r)\"\n# r_nm pg\n";
        s<<"# r_nm rho_O_nm-3 volume_sum_nm3 water_count MC_hits\n";
        int low=0;
        for(size_t k=0;k<counts_.size();++k) {
            double r=(k*bin_+std::min((k+1)*bin_,rmax_))/2;
            double rho=volumes_[k]>0?counts_[k]/volumes_[k]:std::numeric_limits<double>::quiet_NaN();
            o<<r<<" "<<rho/bulkDensity_<<"\n";
            s<<r<<" "<<rho<<" "<<volumes_[k]<<" "<<counts_[k]<<" "<<hits_[k]<<"\n";
            if(hits_[k]<100) ++low;
        }
        o.flush(); s.flush();
        if(!o || !s) GMX_THROW(FileIOError("Failed writing output files."));
        if(low) fprintf(stderr,"WARNING: %d bins have fewer than 100 random hits (including inaccessible bins); inspect the density diagnostics.\n",low);
    }
};
TrajectoryAnalysisModulePointer createPgr() { return TrajectoryAnalysisModulePointer(new Pgr); }
