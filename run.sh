echo "Compiling gsd.c++ for reading GSD file"
g++ -g gsd.c++ -o gsd

echo "Reading file : init_strip.gsd"
./gsd init_strip.gsd > out.txt

echo "Check output in out.txt"
