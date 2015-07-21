#include "myElement.h"


myElement::myElement()
{
  val = 0;
}

myElement::myElement(int _val)
{
  val = _val;
}
    
myElement::~myElement()
{
}

ostream& operator<<(ostream& theStream, const myElement& elem)
{
  theStream << elem.getValue();    
  return theStream;
}

int myElement::operator!=(const myElement &rhs) const
{
  if (val != rhs.getValue())
    {
      return 1;
    }
  return 0;
}

int myElement::operator==(const myElement &rhs) const
{
  if (val == rhs.getValue())
    {
      return 1;
    }
  return 0;
}

myElement myElement::operator=(const myElement &rhs)
{
  val = rhs.getValue();
  return rhs;
}

void myElement::setValue(int _val)
{
  val = _val;
}
