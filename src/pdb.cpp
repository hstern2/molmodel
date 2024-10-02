#include <cctype>
#include "pdb.hpp"
#include "coord.hpp"
#include "out.hpp"

PDBAtomRecord::PDBAtomRecord()
{
  record[0] = name[0] = altLoc[0] = resName[0] = chainID[0] = iCode[0] = '\0';
  resSeq = -1;
}

int PDBAtomRecord::is_empty() const
{
  return strlen(name) == 0;
}

int PDBAtomRecord::is_heterogen() const
{
  return !strcmp(record,"HETATM");
}

void PDBAtomRecord::write(ostream &s, const char *sym, int i, const Cartesian &r) const
{
  char buf[1024];
  snprintf(buf, 256, "%-6s%5d ", strlen(record) > 0 ? record : "ATOM", i);
  if (strlen(name) == 0)
    snprintf(buf+12, 256, " %-3s", sym);
  else
    snprintf(buf+12, 256, "%-4s", name);
  char resSeqStr[5];
  resSeqStr[0] = '\0';
  if (resSeq >= 0)
    snprintf(resSeqStr, 5, "%4d",resSeq);
  snprintf(buf+16, 256, "%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f\n",
	  altLoc,resName,chainID,resSeqStr,iCode,r.x,r.y,r.z);
  s << buf;
}

int PDBAtomRecord::is_alternate_location() const
{
  return strlen(altLoc) > 0 && strcmp(altLoc,"A");
}

int PDBAtomRecord::is_hydrogen() const
{
  insist(strlen(name) > 0);
  if (name[0] == 'H')
    return 1;
  if (strlen(name) > 1 && (isdigit(name[0]) || isspace(name[0])) && name[1] == 'H')
    return 1;
  return 0;
}

static void read_string(const char *buf, char *t, int n)
{
  for (int i = 0; i < n; i++)
    if (!isspace(buf[i]))
      *t++ = buf[i];
  *t = '\0';
}

void PDBAtomRecord::read(const char *buf, Cartesian &r)
{
  read_string(buf, record, 6);
  strncpy(name, buf+12, 4);
  name[4] = '\0';
  read_string(buf+16, altLoc, 1);
  read_string(buf+17, resName, 3);
  read_string(buf+21, chainID, 1);
  sscanf(buf+22, "%4d", &resSeq);
  read_string(buf+26, iCode, 1);
  sscanf(buf+30, "%8lf", &r.x);
  sscanf(buf+38, "%8lf", &r.y);
  sscanf(buf+46, "%8lf", &r.z);
  snprintf(atomID, 32, "%s %s %s %d%s", name, resName, chainID, resSeq, iCode);
}

static void write_str(ostream &s, const char *c)
{
  if (strlen(c) > 0)
    s << c << " ";
  else
    s << ". ";
}

static void read_str(istream &s, char *c)
{
  s >> c;
  if (!strcmp(c,"."))
    *c = 0;
}

istream & operator>>(istream &s, PDBAtomRecord &p)
{
  read_str(s,p.record);
  read_str(s,p.name);
  read_str(s,p.altLoc);
  read_str(s,p.resName);
  read_str(s,p.chainID);
  read_str(s,p.iCode);
  return s >> p.resSeq;
}

ostream & operator<<(ostream &s, const PDBAtomRecord &p)
{
  write_str(s,p.record);
  write_str(s,p.name);
  write_str(s,p.altLoc);
  write_str(s,p.resName);
  write_str(s,p.chainID);
  write_str(s,p.iCode);
  return s << p.resSeq;
}

Str PDBAtomRecord::symbol() const
{
  char buf[3];
  if (isdigit(name[0]) || isspace(name[0])) {
    buf[0] = name[1];
    buf[1] = '\0';
  } else {
    buf[0] = name[0];
    buf[1] = tolower(name[1]);
  }
  buf[2] = '\0';
  return Str(buf);
}
