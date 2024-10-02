#ifndef TIMING_H
#define TIMING_H

#ifdef TIMING
#define TIMESTART(x) tmstart(x)
#define TIMESTOP(x) tmstop(x)
#define TIMEWRITE() tmwrite()
#else
#define TIMESTART(x) ((void) 0)
#define TIMESTOP(x) ((void) 0)
#define TIMEWRITE() ((void) 0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

  void tmstart(const char *s);
  void tmstop(const char *s);
  void tmwrite();

#ifdef __cplusplus
}
#endif

#endif /* TIMING_H */
