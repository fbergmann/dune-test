// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_COPASI_CONTROLLER_HH
#define DUNE_COPASI_CONTROLLER_HH

/** \brief Control time step after reaction.

    If some concentration is negative, then it returns false.
    Otherwise it returns true (works in parallel)

    \tparam V           Vector backend (containing DOF)
*/
template<class GV, class V>
bool controlReactionTimeStepParallel (GV& gv, V& v)
{
  int passed = 1;
  for (auto it=v.begin();it!=v.end();++it)
    if (*it < 0.)
      {
        passed = 0;
        break;
      }
  if (gv.comm().size()>1)
    passed =  gv.comm().min(passed);
  if (passed) return true;
  else return false;
}

#endif
