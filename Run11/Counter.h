#ifndef COUNTER_H
#define COUNTER_H

struct Counter
{
  Counter(){
  }
  Counter(int totNum, int totAccept){
	  totalNumber = totNum;
	  totalAccepted = totAccept;
  }
  int totalNumber = 0;
  int totalAccepted = 0;
  
  double getRatio() const
  {
    return static_cast<double>(totalAccepted)/totalNumber;
  }

  //  std::cout<<"totalNumber: "<<totalNumber<<std::endl;
  // std::cout<<"totalAccepted: "<<totalAccepted<<std::endl;
  
};
#endif /* !COUNTER_H */
