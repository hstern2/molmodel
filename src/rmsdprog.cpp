#include "cvec.hpp"
#include "rmsd.hpp"

const char *usage =
"Usage: rmsd <file1> <file2>\n"
"\n"
"Aligns the coordinates in <file2> with those in <file1>.\n"
"Prints out the root-mean-square distance, and the aligned coordinates.\n"
"<file1> and <file2> should contain a list of coordinates\n"
"separated by whitespace, with the same number in each\n";

int main(int argc, char *argv[])
{
  if (argc != 3) {
    printf("%s", usage);
    exit(1);
  }
  Lst<Cartesian> v1, v2;
  istream *f1 = StreamSearch(argv[1]), *f2 = StreamSearch(argv[2]);
  Cartesian tmp;
  while (*f1 >> tmp)
    v1.add(tmp);
  while (*f2 >> tmp)
    v2.add(tmp);
  CVec c1, c2;
  c1.copy_from_list(v1);
  c2.copy_from_list(v2);
  if (c1.size() != c2.size())
    die("rmsd: files must contain the same number of coordinates");
  printf("RMSD: %f\n", RootMeanSquareDistance(c1,c2));
  printf("Second set of coordinates, aligned to first set:\n");
  for (int i = 0; i < v2.size(); i++)
    printf("%15.8f %15.8f %15.8f\n", c2[i].x, c2[i].y, c2[i].z);
}

