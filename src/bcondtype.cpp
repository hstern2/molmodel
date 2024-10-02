#include <cstdlib>
#include <cstring>
#include "bcondtype.h"
#include "out.hpp"
#include "str.hpp"

istream & operator>>(istream &stream_ref, bcondtype_t &enum_obj)
{
  Str buf;
  stream_ref >> buf;
  if (!strcmp(buf, "non_periodic"))
    enum_obj = non_periodic;
  else if (!strcmp(buf, "cubic"))
    enum_obj = cubic;
  else if (!strcmp(buf, "orthorhombic"))
    enum_obj = orthorhombic;
  else if (!strcmp(buf, "bcc"))
    enum_obj = bcc;
  else if (!strcmp(buf, "fcc"))
    enum_obj = fcc;
  else if (!strcmp(buf, "triclinic"))
    enum_obj = triclinic;
  else
    die("bcondtype_t: must be one of 'non_periodic','cubic','orthorhombic','bcc','fcc','triclinic'");
  return stream_ref;
}

ostream & operator<<(ostream &stream_ref, const bcondtype_t &enum_obj)
{
  switch (enum_obj) {
  case non_periodic:
    return stream_ref << "non_periodic";
  case cubic:
    return stream_ref << "cubic";
  case orthorhombic:
    return stream_ref << "orthorhombic";
  case bcc:
    return stream_ref << "bcc";
  case fcc:
    return stream_ref << "fcc";
  case triclinic:
    return stream_ref << "triclinic";
  }
  return stream_ref;
}
