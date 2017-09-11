ulimit -c unlimited


echo "Compiling gsd.c++ for reading GSD file"
g++ -g gsd_energy.c++ -o gsd

echo "Reading file : init_strip.gsd"
./gsd ../Sim_dump_wall/init_strip.gsd $1 > ../Sim_dump_wall/gsd_out.txt

echo "Check output in ../Sim_dump/gsd_out.txt"
