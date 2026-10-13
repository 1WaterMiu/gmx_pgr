# Source this file in Bash: source ./activate.sh
_pgr_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${_pgr_root}/bin:${PATH}"
unset _pgr_root
hash -r
