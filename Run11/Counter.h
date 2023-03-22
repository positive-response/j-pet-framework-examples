#ifndef COUNTER_H
#define COUNTER_H

struct Counter
{
  Counter(){
  }
  int totalNumber = 0;
  int totalAccepted = 0;
  
  double getRatio() const
  {
    return static_cast<double>(totalAccepted)/totalNumber;
  }
  
};
#endif /* !COUNTER_H */
