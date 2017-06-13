#ifndef CUTS_ENUM_H_
#define CUTS_ENUM_H_

#include <string>
#include <functional>
template< typename T >
class Enum {
public:
  class Iterator {

  public:
    Iterator( int value ) :  m_value( value ) { }

    T operator*( void ) const {
      return (T)m_value;
    }

    void operator++( void ) {
      ++m_value;
    }

    bool operator!=( Iterator rhs ) {
      return m_value != rhs.m_value;
    }

  private:
    int m_value;
  };

};

template< typename T >
typename Enum<T>::Iterator begin( Enum<T> ) {
  return typename Enum<T>::Iterator( (int)T::First );
}

template< typename T >
typename Enum<T>::Iterator end( Enum<T> ) {
  return typename Enum<T>::Iterator( ((int)T::Last) + 1 );
}


struct EnumHash {
  template<typename T> inline typename std::enable_if<std::is_enum<T>::value, std::size_t>::type

  operator()(const T&t) const  {
    return static_cast<std::size_t>(t);
  }

};


enum class CUTS {
  eGen,
  eGTau,        eGTop,        eGElec,       eGMuon,       eGZ,        eGW,       eGHiggs,
  eRVertex,     eRMuon1,      eRMuon2,      eRElec1,      eRElec2,    eRTau1,   eRTau2,
  eRJet1,       eRJet2,       eRCenJet,     eR1stJet,     eR2ndJet,   eRBJet,   eRWjet,
  eDiElec,      eDiMuon,      eDiTau,       eDiJet,
  eMuon1Tau1,   eMuon1Tau2,   eMuon2Tau1,   eMuon2Tau2,
  eElec1Tau1,   eElec1Tau2,   eElec2Tau1,   eElec2Tau2,
  eMuon1Elec1,  eMuon1Elec2,  eMuon2Elec1,  eMuon2Elec2,
  eSusyCom,     eMET,         eNuTau,       eRTrig1,      eRTrig2,
  First = eGen,
  Last = eRTrig2};


  #endif
