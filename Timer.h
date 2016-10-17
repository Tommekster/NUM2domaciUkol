#pragma once


class Timer
{
   public:
   
      Timer();

      void reset();

      void stop();

      void start();

      double getRealTime();

      double getCPUTime();

      unsigned long long int getCPUCycles();      
         
   protected:

   double initialRealTime, totalRealTime,
          initialCPUTime, totalCPUTime;
   
   unsigned long long int initialCPUCycles, totalCPUCycles;
   
   bool stopState;
   
   inline unsigned long long rdtsc()
   {
     unsigned hi, lo;
     __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
     return ( ( unsigned long long ) lo ) | ( ( ( unsigned long long ) hi ) << 32 );
   }
};

