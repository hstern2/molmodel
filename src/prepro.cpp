#include <cstring>
#include <cctype>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "prepro.hpp"
#include "out.hpp"
#include "lst.hpp"

PreprocessedStream::PreprocessedStream() :
  ftmp(0), findex(0), s(0)
#ifdef USE_MPI
  , size(0), rank(0)
#endif
{ }

void PreprocessedStream::init(int n, char **fnames)
{
  cleanup();
  n--;
  fnames++;
  /* Set variable 'pwd' to basename of current working directory */
  char tmpbuf[1024];
  getcwd(tmpbuf, 1024);
  char *b;
  for (b = &tmpbuf[strlen(tmpbuf) - 1]; b >= tmpbuf; b--)
    if (*b == '/')
      break;
  var["pwd"] = ++b;
  /* Some more predefined variables:  the process ID and MPI rank */
  const int pid = getpid();
  snprintf(tmpbuf,1024,"%d",pid);
  var["pid"] = tmpbuf;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Str leading_zero = "";
  if (size >= 100) {
    if (rank < 10)
      leading_zero = "00";
    else if (rank < 100)
      leading_zero = "0";
  } else if (size >= 10 && rank < 10) {
    if (rank < 10)
      leading_zero = "0";
  }
  char pbuf[256];
  sprintf(pbuf,"%s%d",(const char *) leading_zero, rank);
  var["rank"] = pbuf;
#endif
  /* Create a temporary and index file */
  TempName(tmpfname);
  ftmp = FileStream(tmpfname);
  TempName(index_fname);
  findex = FileStream(index_fname);
  bool anyfile = false;
  /* Loop through command line arguments */
  while (n > 0)
    if (strlen(*fnames) > 1 && (*fnames)[0] == '-') {
      /* Variable definition */
      const char *vname = &(*fnames)[1];
      fnames++, n--;
      if (n == 0)
	die("PreprocessedStream: expecting variable definition for '-%s'", vname);
      var[vname] = *fnames;
      fnames++, n--;
    } else {
      /* Input file */
      anyfile = true;
      Str tmp = *fnames;
      for (char *c = ((char *) tmp) + strlen(tmp) - 1; c > (char *) tmp; c--)
	if (*c == '.') {
	  *c = '\0';
	  break;
	} 
      var["base"] = tmp;
      read_fname(*fnames);
      fnames++, n--;
    }
  if (!anyfile)
    read_file(stdin,"stdin");
  if (ftmp)
    delete ftmp;
  ftmp = 0;
  if (findex)
    delete findex;
  findex = 0;
  s = StreamSearch(tmpfname);
  if (!s)
    die("PreprocessedStream: could not open temporary file");
}

PreprocessedStream::~PreprocessedStream()
{ cleanup(); }

void PreprocessedStream::cleanup()
{
  if (strlen(tmpfname) > 0)
    remove(tmpfname);
  if (strlen(index_fname) > 0)
    remove(index_fname);
  tmpfname = index_fname = "";
  if (s)
    delete s;
  s = 0;
}

void PreprocessedStream::read_fname(const char *fname)
{
  FILE *f = FileSearch(fname);
  read_file(f,fname);
  fclose(f);
}

void PreprocessedStream::read_file(FILE *f, const char *fname, bool expecting_end)
{
  char ln[1024];
  int nline = 0;
  bool ignore = false;
  while(fgets(ln, 1024, f)) {
    nline++;
    if (strlen(ln) == 1024)
      die("%s:%d: PreprocessedStream::read_file: line too long", fname, nline);
    char *c;
    for (c = ln; *c; c++)
      if (*c == '#') { // comment
	*c++ = '\n';
	*c = '\0';
      }
    for (c = ln; isspace(*c); c++)
      if (*c == '!')
	break;
    if (!strncmp(c, "!ignore", 7))
      ignore = true;
    if (ignore) {
      if (!strncmp(c, "!end", 4)) {
	ignore = false;
	for (char *ctmp = ln; ctmp < c+4; ctmp++)
	  *ctmp = ' ';
      } else {
	continue;
      }
    }
    if (*c != '!' || 
	!strncmp(c, "!define", 7) || 
	!strncmp(c, "!parallel_define", 16) ||
	!strncmp(c, "!include", 8) ||
	!strncmp(c, "!path", 5)) {
      read_line(ln,fname,nline); // ordinary line, define, include
      continue;
    } else if (!strncmp(c, "!end", 4)) {
      if (!expecting_end)
	die("%s:%d: PreprocessedStream::read_line: unexpected !end", fname, nline);
      return;
    }
    if (!strncmp(c, "!foreach", 8)) {
      istringstream s(c + 8);
      Str v;
      Lst<Str> vals;
      s >> v >> vals;
      if (!s || vals.size() == 0)
	die("%s:%d: PreprocessedStream::expecting !foreach var { values ... }", fname, nline);
      LstItr<Str> vv;
      fpos_t pos;
      fgetpos(f,&pos);
      for (vv.init(vals); vv.ok(); vv.next()) {
	var[v] = vv();
	fsetpos(f,&pos);
	read_file(f,fname,1);
      }
      var[v] = "";
    } else if (!strncmp(c, "!repeat", 7)) {
      istringstream s(c + 7);
      int n;
      s >> n;
      if (!s || n <= 0)
	die("%s:%d: PreprocessedStream::expecting !repeat n where n > 0",fname,nline);
      fpos_t pos;
      fgetpos(f,&pos);
      for (int i = 0; i < n; i++) {
	fsetpos(f,&pos);
	read_file(f,fname,1);
      }
    } else {
      die("%s:%d: PreprocessedStream::read_file: unknown command: %s",fname,nline,ln);
    }
  }
  if (expecting_end)
    die("%s:%d: PreprocessedStream::read_file: expecting !end, reached end of file",fname,nline);
}

void PreprocessedStream::read_line(char *ln, const char *fname, int nline) 
{
  char *c, *c2;
  /* Look for variable replacements */
  for (c = ln; *c; c++)
    if (*c == '$') {
      *c++ = '\0';
      Str tmp(ln);
      char vname[256];
      int i;
      for (i = 0; (isalnum(*c) || *c == '_') && i < 256; i++, c++)
	vname[i] = *c;
      if (i >= 256)
	die("%s:%d: PreprocessedStream::variable name too long",fname,nline);
      vname[i] = '\0';
      if (strlen(vname) == 0)
	die("%s:%d: PreprocessedStream::read_file: no variable name after '$'",fname,nline);
      if (var.exists(vname) && strlen(var[vname]) > 0)
	tmp += var[vname];
      else if ((c2 = getenv(vname)) != 0)
	tmp += c2;
      else
	die("%s:%d: PreprocessedStream::read_file: undefined variable: '%s'",fname,nline,vname);
      tmp += c;
      read_line(tmp,fname,nline); // process the modified line again
      return;
    }
  /* Look for path */
  for (c = ln; *c; c++)
    if (!strncmp(c, "!path", 5)) {
      *c++ = '\n';
      *c++ = '\0';
      read_line(ln,fname,nline); // process line before path directive
      c += 3;
      AddToSearchPath(c);
      return;
    }
  /* Look for variable definitions */
  for (c = ln; *c; c++)
    if (!strncmp(c, "!define", 7)) {
      *c++ = '\n';
      *c++ = '\0';
      read_line(ln,fname,nline); // process the line before the define
      c += 5;
      char *vname = strtok(c, " \t\n");
      if (!vname)
	die("%s:%d: PreprocessedStream::read_file: no variable name after !define",fname,nline);
      for (char *tmp = vname; *tmp; tmp++)
	if (!(isalnum(*tmp) || *tmp == '_'))
	  die("%s:%d: PreprocessedStream::read_file: illegal variable name: %s",fname,nline,vname);
      while(*c++)
	;
      for (char *d = c + strlen(c) - 1; d > c && isspace(*d); d--)
	*d = '\0'; // remove trailing whitespace
      var[vname] = c;
      return;
    }
  /* Look for parallel variable definitions */
  for (c = ln; *c; c++)
    if (!strncmp(c, "!parallel_define", 16)) {
#ifdef USE_MPI
      *c++ = '\n';
      *c++ = '\0';
      read_line(ln,fname,nline); // process the line before the define
      c += 14;
      istringstream s(c);
      Str vname;
      s >> vname;
      if (!s)
	die("%s:%d: PreprocessedStream::parallel_define: expecting variable name", fname, nline);
      Lst<Str> vals;
      s >> vals;
      if (!s || vals.size() != size)
	die("%s:%d: PreprocessedStream::parallel_define: expecting list of %d "
	    "values (that's how many processes we are running on", fname, nline, size);
      LstItr<Str> vi;
      int i;
      for (i = 0, vi.init(vals); vi.ok(); i++, vi.next())
	if (i == rank) {
	  var[vname] = vi();
	  return;
	}
#else
      die("%s:%d PreprocessedStream: this executable was not compiled with MPI "
	  "so parallel_define is not available",fname,nline);
#endif
    }
  /* Look for includes */
  int quote = 0;
  for (c = ln; *c; c++) {
    if (*c == '"')
      quote = !quote;
    if (quote)
      continue;
    if (!strncmp(c, "!include ", 9) || *c == '<') {
      if (*c == '<') 
	*c++ = '\0';
      else {
	*c++ = '\0';
	c += 8;
      }
      while (*c && isspace(*c))
	c++;
      read_line(Str(ln),fname,nline);
      read_line(Str(" "),fname,nline);
      char *d = new char[strlen(c)+1];
      int nc = 0;
      while (*c && !isspace(*c) && *c != '{' && *c != '}')
	d[nc++] = *c++;
      d[nc] = '\0';
      if (strlen(d) == 0)
	die("%s:%d: PreprocessedStream::read_file: expecting filename after '!include' or '<'",fname,nline);	
      read_fname(d);
      read_line(c,fname,nline);
      delete[] d;
      return;
    }
  }
  /* Nothing to do to this line -- write it to the temp file */
  if (!(*ftmp << ln))
    die("PreprocessedStream::read_line: cannot write to temp file (disk full?)");
  /* Write index file */
  if (strcmp(last_fname, fname) || nline != last_nline) {
    *findex << ftmp->tellp() << " " << fname << " " << nline << "\n";
    if (strcmp(last_fname, fname))
      last_fname = fname;
    last_nline = nline;
  }
}

void PreprocessedStream::show_self()
{
  FILE *f = fopen(tmpfname, "r");
  if (!f) {
    Out() << "PreprocessedStream::show_self: cannot open temporary file.\n";
    return;
  }
  int c;
  while ((c = getc(f)) != EOF)
    OutPrintf("%c",c);
  Out() << flush;
  fclose(f);
}

void PreprocessedStream::write_current_file_and_line(ostream &stmp) const
{
  FILE *f = fopen(index_fname,"r");
  if (!(s && f))
    return;
  const int pos = s->tellg();
  char buf[1024], fname[1024];
  int mypos, nline;
  while (fgets(buf, 1024, f)) {
    sscanf(buf,"%d%s%d",&mypos,fname,&nline);
    if (mypos >= pos) {
      stmp << fname << ":" << nline << "\n" << flush;
      break;
    }
  }
  fclose(f);
}
