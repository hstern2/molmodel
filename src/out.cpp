#include <fstream>
#include <unistd.h>
#include <cstdlib>
#include <dirent.h>
#include <sys/types.h>
#include <cstdarg>
#include "out.hpp"
#include "lst.hpp"

static ostream *out = 0;
static int precision = 8;
static Lst<Str> *search_path = 0;

void AddToSearchPath(const char *p)
{
  if (p) {
    if (!search_path)
      search_path = new Lst<Str>;
    search_path->add(p);
  }
}

static bool is_directory(const char *fname)
{
  DIR *dir = opendir(fname);
  if (dir) {
    closedir(dir);
    return true;
  }
  return false;
}

FILE *FileSearch(const char *fname)
{
  if (strlen(fname) == 0)
    die("FileSearch: expecting a filename");
  FILE *f = fopen(fname, "r");
  if (f && !is_directory(fname))
    return f;
  Str fn(fname);
  if (search_path) {
    LstItr<Str> s;
    for (s.init(*search_path); s.ok(); s.next()) {
      const Str tmpf = s() + "/" + fn;
      f = fopen(tmpf,"r");
      if (f && !is_directory(tmpf))
	return f;
    }
  }
  die("FileSearch: cannot open %s", (const char *) fn);
  return 0;
}

istream *StreamSearch(const char *fname)
{
  istream *f = new ifstream(fname);
  if (f && *f && !is_directory(fname))
    return f;
  Str fn(fname);
  if (search_path) {
    LstItr<Str> s;
    for (s.init(*search_path); s.ok(); s.next()) {
      Str d = s() + "/" + fn;
      f = new ifstream((const char *) d);
      if (f && *f && !is_directory(d))
	return f;
    }
  }
  die("FileSearch: cannot open %s", (const char *) fn);
  return 0;
}

void TempName(Str &s)
{
  char tmp[256];
  int n = getpid();
  static int noffset = 0;
  FILE *f;
  do
    snprintf(tmp, 256, "/tmp/tmp%d.%d", n, noffset++);
  while ((f = fopen(tmp, "r")) && !fclose(f));
  s = tmp;
}

void OutInit() 
{ 
  SetPrecision(precision);
}

void OutFinalize()
{
  if (search_path) {
    delete search_path;
    search_path = 0;
  }
}

void OutSetToFile(const char *f) 
{
  if (out)
    delete out;
  out = FileStream(f);
}

void SetPrecision(int p)
{
  Out().precision(precision = p);
}

ostream &Out() { return out ? *out : cout; }

ostream *FileStream(const char *fname)
{
  ostream *f = new ofstream(fname);
  if (!(f && *f)) {
    Out() << "FileStream: cannot create " << fname << "\n" << flush;
    die("");
  }
  f->precision(precision);
  return f;
}

ostream *AppendFileStream(const char *fname)
{
  ostream *f = new ofstream(fname, ios_base::app);
  if (!(f && *f)) {
    Out() << "FileStream: cannot append to " << fname << "\n" << flush;
    die("");
  }
  f->precision(precision);
  return f;
}

bool end_of_file(FILE *f)
{
  int c = fgetc(f);
  if (c == EOF)
    return true;
  ungetc(c, f);
  return false;
}

static char spaces[21] = "                    ";
static int spi = 20;

void IndentPush()
{
  spi -= 2;
}

void IndentPop()
{
  spi += 2;
}

const char *Indent()
{
  int i = spi;
  if (i > 20)
    i = 20;
  else if (i < 0)
    i = 0;
  return spaces + i;
}

void OutPrintf(const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  char *ret;
  vasprintf(&ret, fmt, argp);
  va_end(argp);
  Out() << ret;
  free(ret);
}

void die(const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  char *ret;
  vasprintf(&ret, fmt, argp);
  va_end(argp);
  Out() << ret << "\n" << flush;
  free(ret);
  throw FatalException();
}
