#include <unistd.h>
#include <ctime>
#include "prog.hpp"
#include "lst.hpp"
#include "str.hpp"
#include "timing.h"
#include "out.hpp"

void Program::exit()
{
  die("Called Program::exit.");
}

void Program::system(istream &s)
{
  char c[1024];
  int i = 0;
  while (s.get(c[i]) && c[i] != '\n' && i < 1024)
    i++;
  c[i] = '\0';
  ::system(c);
}

void Program::write_header()
{
  /* Get the date/time and hostname */
  time_t t = time(0);
  char buf[256];
  gethostname(buf, 256);
  Out() << ctime(&t) << buf << "\n"
	<< "Process id: " << getpid() << "\n\n";
#ifdef SHOW_RCSID
  << "To check out RCS source code files for this executable:\n";
  for (int i = 0; i < nRcsId; i++)
    Out() << RcsId[i];
  Out() << "\n" << flush;
#endif
}

void Program::timing()
{
  TIMEWRITE();
}

void Program::output_precision(int n)
{
  if (n < 0)
    die("Program::output_precision: expecting a positive int");
  SetPrecision(n);
}

void Program::output_file(Str f)
{
  if (!strcmp(f, "stdout"))
    OutInit();
  else
    OutSetToFile(f);
}

void Program::add_to_search_path(Str f)
{
  AddToSearchPath(f);
}
