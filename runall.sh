nuc(){
  python3 init.py "./dataset/$1"
}
for f in $(ls ./dataset); do nuc "$f" & done
