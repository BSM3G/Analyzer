#ifndef _FILL_INFO_H_
#define _FILL_INFO_H_

enum class FILLER { Single, Dipart, Dilepjet, None };

struct FillVals {
  CUTS ePos;
  FILLER type;
  Particle* part;
  Particle* part2;

  FillVals(): type(FILLER::None) {}
  FillVals(CUTS _ePos): ePos(_ePos), type(FILLER::None), part(nullptr), part2(nullptr) {}
  FillVals(CUTS _ePos, FILLER _type, Particle* _part, Particle* _part2=nullptr) :
  ePos(_ePos), type(_type), part(_part), part2(_part2) { }

};

struct GenFill {
  int status;
  CUTS ePos;
  GenFill(int _status, CUTS _ePos) : status(_status), ePos(_ePos) {}
};


#endif
