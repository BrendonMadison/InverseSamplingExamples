#include <vector>

using namespace std;

//Returns a std vector of gaussian random variables using the analytical solution
std::vector<Double_t> Gaus(Double_t MU, Double_t VAR, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> gRan;

  for(Int_t i = 0; i < NUM; i++)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      gRan.push_back(TMath::ErfInverse(rg->Uniform() * 2.0 - 1.0) * sqrt(2.0)*VAR + MU);
    }
  delete rg;
  delete ts;
  return gRan;
}

//Returns a std vector of gaussian random variable but times it for comparison sake
int GausTime()
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  TTimeStamp *cs = new TTimeStamp();
  Double_t gRan = 0.0;
  cs->Set();
  Double_t initime = cs->AsDouble();
  for(Int_t i = 0; i < 1000000; i++)
    {
      //      ts->Set();
      //      rg->SetSeed(ts->GetNanoSec()); Don't need this if we are comparing to method that doesn't have dynamic seeding
      gRan = TMath::ErfInverse(rg->Uniform() * 2.0 - 1.0) * sqrt(2.0)*1.0;
    }
  cs->Set();
  cout << "Time to generate 1M Gaussian Variates: " << cs->AsDouble() - initime << endl;
  delete rg;
  delete ts;
  delete cs;
  return 0;
}

//Returns a std vector of exponential random variables using the analytical solution
std::vector<Double_t> Expo(Double_t MU, Double_t LAM, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> eRan;

  for(Int_t i = 0; i < NUM; i++)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      eRan.push_back(-1.0/LAM * log(1.0 - rg->Uniform()) + MU);
    }
  delete rg;
  delete ts;
  return eRan;
}

//Returns a std vector of exponential distribution of a gaussian random variable using the analytical solutions
//This is also known as the exponentially modified gaussian, but I use a different notation from typical 
//https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
//I verified that it is as I checked the mode, variance values and they match
std::vector<Double_t> ExpoGaus(Double_t MU, Double_t LAM, Double_t VAR, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> egRan;
  Double_t tmp;

  Int_t cnt = 0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = TMath::ErfInverse(rg->Uniform() * 2.0 - 1.0) * sqrt(2.0)*VAR + MU;
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp += -1.0/LAM * log(1.0 - rg->Uniform());
      if (TMath::IsNaN(tmp) == false)
	{
	  egRan.push_back(tmp - VAR/LAM);
	  cnt++;
	}
    }
  delete rg;
  delete ts;
  return egRan;
}

//Returns a Cauchy distribution random variate
std::vector<Double_t> Cauchy(Double_t MU, Double_t GAM, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = GAM*(tan(3.1415926*(rg->Uniform() - 0.5) ) + MU);
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	}
    }
  delete rg;
  delete ts;
  return xRan;
}

//Returns a Levy distribution random variate
std::vector<Double_t> Levy(Double_t MU, Double_t KAP, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = MU + KAP/(2.0*pow(TMath::ErfInverse(1.0 - rg->Uniform()),2.0));
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	}
    }
  delete rg;
  delete ts;
  return xRan;
}

//Returns a arcsine distribution random variate on the bounds of  A <= x <= B
std::vector<Double_t> ArcSin(Double_t A, Double_t B, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = pow(sin(3.1415926 * rg->Uniform() * 0.5),2.0) * (B - A) + A;
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	}
    }
  delete rg;
  delete ts;
  return xRan;
}

//Returns a skewed arcsine distribution random variate on the bounds of  A <= x <= B
//Uses the M parameter, which is the power of the skew... not intuitive!
std::vector<Double_t> SkewArcSin(Double_t A, Double_t B, Double_t M, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  Double_t Lcnt = 0.0;
  Double_t Rcnt = 0.0;
  Double_t Alr = 0.0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = pow(sin(3.1415926 * rg->Uniform() * 0.5),1.0/M) * (B - A) + A;
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	  if (tmp < (B - A)/2.0)
	    {
	      Lcnt += 1.0;
	    }
	  else
	    {
	      Rcnt += 1.0;
	    }
	}
    }
  Alr = (Lcnt - Rcnt)/(Lcnt + Rcnt);
  cout << "Left-right asymmetry: " << Alr << endl;
  delete rg;
  delete ts;
  return xRan;
}

//Returns a skewed arcsine distribution random variate on the bounds of  A <= x <= B
//This is merely for checking the Left-Right Asymmetry of the skewed arcsine (A_LR)
Double_t GetAlr(Double_t A, Double_t B, Double_t M, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  Double_t tmp;

  Int_t cnt = 0;
  Double_t Lcnt = 0.0;
  Double_t Rcnt = 0.0;
  Double_t Alr = 0.0;
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = pow(sin(3.1415926 * rg->Uniform() * 0.5),1.0/M) * (B - A) + A;
      if (TMath::IsNaN(tmp) == false)
	{
	  cnt++;
	  if (tmp < (B - A)/2.0)
	    {
	      Lcnt += 1.0;
	    }
	  else
	    {
	      Rcnt += 1.0;
	    }
	}
    }
  Alr = (Lcnt - Rcnt)/(Lcnt + Rcnt);
  //cout << "Left-right asymmetry: " << Alr << endl;
  delete rg;
  delete ts;
  return Alr;
}

//Returns a skewed arcsine distribution random variate on the bounds of  A <= x <= B using A_LR
//This one uses A_LR (as PLR) which is far more intuitive!
std::vector<Double_t> ArcSinALR(Double_t A, Double_t B, Double_t PLR, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  Double_t Lcnt = 0.0;
  Double_t Rcnt = 0.0;
  Double_t Alr = 0.0;
  Double_t PI = 3.141592653589793;
  Double_t M1 = log((PLR + 1.0)/(1.40529))*-1.0/0.747437;
  Double_t M2 = (1.0/3.22864) * ((1.44148/(PLR + 0.561397)) - 1.0);
  Double_t M = 0.0;
  if (M1 > 0.6)
    {
      M = M1;
    }
  else
    {
      M = M2;
    }
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = pow(sin(PI * rg->Uniform() * 0.5),1.0/M) * (B - A) + A;
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	  if (tmp < (B - A)/2.0)
	    {
	      Lcnt += 1.0;
	    }
	  else
	    {
	      Rcnt += 1.0;
	    }
	}
    }
  Alr = (Lcnt - Rcnt)/(Lcnt + Rcnt);
  cout << "Left-right asymmetry: " << Alr << endl;
  delete rg;
  delete ts;
  return xRan;
}

//Returns a skewed arcsine distribution random variate on the bounds of  A <= x <= B using A_LR
//Used for sanity checking!
Double_t CheckALR(Double_t A, Double_t B, Double_t PLR, Int_t NUM)
{
  TRandom3 *rg = new TRandom3();
  TTimeStamp *ts = new TTimeStamp();
  std::vector<Double_t> xRan;
  Double_t tmp;

  Int_t cnt = 0;
  Double_t Lcnt = 0.0;
  Double_t Rcnt = 0.0;
  Double_t Alr = 0.0;
  Double_t PI = 3.141592653589793;
  Double_t M1 = log((PLR + 1.0)/(1.40529))*-1.0/0.747437;
  Double_t M2 = (1.0/3.22864) * ((1.44148/(PLR + 0.561397)) - 1.0);
  Double_t M = 0.0;
  if (M1 > 0.6)
    {
      M = M1;
    }
  else
    {
      M = M2;
    }
  while(cnt < NUM)
    {
      ts->Set();
      rg->SetSeed(ts->GetNanoSec());
      tmp = pow(sin(PI * rg->Uniform() * 0.5),1.0/M) * (B - A) + A;
      if (TMath::IsNaN(tmp) == false)
	{
	  xRan.push_back(tmp);
	  cnt++;
	  if (tmp < (B - A)/2.0)
	    {
	      Lcnt += 1.0;
	    }
	  else
	    {
	      Rcnt += 1.0;
	    }
	}
    }
  Alr = (Lcnt - Rcnt)/(Lcnt + Rcnt);
  cout << "Left-right asymmetry: " << Alr << endl;
  delete rg;
  delete ts;
  return Alr;
}
