/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "StackOperators.h"
#include "Stackwavefunction.h"
#include "Stackspinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "pario.h"
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>


long SpinAdapted::getRequiredMemoryForWavefunction(const StateInfo& sr, const StateInfo& sc, const std::vector<SpinQuantum>& q) {

  long memory = 0;

  for (int i=0; i < sr.quanta.size(); i++)
    for (int j=0; j<sc.quanta.size(); j++) {
      bool allowedcoupling = false;
      for (int k=0; k<q.size(); k++) {
	if (q[k].allow(sr.quanta[i], sc.quanta[j])) {
	  allowedcoupling = true;
	  break;
	}
      }
      if (allowedcoupling)
	memory += sr.quantaStates[i]*sc.quantaStates[j]+ CACHEBUFFER;
    }

  return memory;
}

void SpinAdapted::StackWavefunction::copyData(const StackWavefunction& a) {
  for (int i=0; i<a.nonZeroBlocks.size(); i++) {
    int lQ = a.nonZeroBlocks[i].first.first, rQ = a.nonZeroBlocks[i].first.second;
    //assert(lQ == nonZeroBlocks[i].first.first && rQ == nonZeroBlocks[i].first.second);
    
    copy(a.nonZeroBlocks[i].second, operator()(lQ, rQ));
   }
   return ;

 }


 double* SpinAdapted::StackWavefunction::allocateWfnOperatorMatrix()
 {
   long index = 0;
   for (int i=0; i<nonZeroBlocks.size(); i++) {
     int lQ = nonZeroBlocks[i].first.first, rQ = nonZeroBlocks[i].first.second;
     StackMatrix m(&data[index], operator()(lQ, rQ).Nrows(), operator()(lQ, rQ).Ncols());
     //operatorMatrix(lQ,rQ).allocate(&data[index], operatorMatrix(lQ, rQ).Nrows(), operatorMatrix(lQ, rQ).Ncols());
     index += operator()(lQ, rQ).Nrows()* operator()(lQ, rQ).Ncols()+ CACHEBUFFER;
     nonZeroBlocks[i].second= m;
     mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), i));
   }
   totalMemory = index;
   return &data[index];
 }

 void SpinAdapted::StackWavefunction::initialise(const StackWavefunction& w)
 {
   *this = w;
   data = Stackmem[omprank].allocate(totalMemory);

   long index = 0;
   for (int i = 0; i<nonZeroBlocks.size(); i++) {
     int lQ = nonZeroBlocks[i].first.first, rQ = nonZeroBlocks[i].first.second;
     StackMatrix m(&data[index], operator()(lQ, rQ).Nrows(), operator()(lQ, rQ).Ncols());
     //operatorMatrix(lQ,rQ).allocate(&data[index], operatorMatrix(lQ, rQ).Nrows(), operatorMatrix(lQ, rQ).Ncols());
     index += operator()(lQ, rQ).Nrows()* operator()(lQ, rQ).Ncols()+ CACHEBUFFER;
     nonZeroBlocks[i].second= m;
     mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), i));
   }
   if (index != totalMemory) exit(0);
 }

 void SpinAdapted::StackWavefunction::initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_)
 {
   long totalMemory = getRequiredMemoryForWavefunction(sl, sr, dQ);
   double* data = Stackmem[omprank].allocate(totalMemory);
   //memset(data, 0, totalMemory*sizeof(double));
   initialise(dQ, sl, sr, onedot_, data, totalMemory);  
 }


 void SpinAdapted::StackWavefunction::initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_, double* pData, long ptotalMemory)
 {
   // initialized a ket wavefunction
   //operatorMatrix.resize(0,0);
   data = pData;
   totalMemory = ptotalMemory;
   initialised = true;
   fermion = false;
   built = true;
   deltaQuantum = dQ;
   onedot = onedot_; 


   nonZeroBlocks.resize(0);
   mapToNonZeroBlocks.clear();
   rowCompressedForm.clear();rowCompressedForm.resize(sl.quanta.size(), vector<int>());
   colCompressedForm.clear();colCompressedForm.resize(sr.quanta.size(), vector<int>());
   //operatorMatrix.resize(sl.quanta.size (), sr.quanta.size ());
   allowedQuantaMatrix.resize(sl.quanta.size (), sr.quanta.size ());

   long index = 0;
   for (int lQ = 0; lQ < sl.quanta.size (); ++lQ)
     for (int rQ = 0; rQ < sr.quanta.size (); ++rQ) {
       allowedQuantaMatrix(lQ, rQ) = false;
       for (int i = 0; i < deltaQuantum.size(); ++i)
	 if (deltaQuantum[i].allow(sl.quanta [lQ] , sr.quanta [rQ])) {
	   allowedQuantaMatrix(lQ, rQ) = true;
	   rowCompressedForm[lQ].push_back(rQ);
	   colCompressedForm[rQ].push_back(lQ);
	   break;
	 }

       if (allowedQuantaMatrix(lQ, rQ))
       {
	 StackMatrix m(&data[index], sl.quantaStates [lQ], sr.quantaStates [rQ]);
	 //operatorMatrix(lQ,rQ).allocate(&data[index], sl.quantaStates [lQ], sr.quantaStates [rQ]);
	 index += sl.quantaStates [lQ]* sr.quantaStates [rQ]+ CACHEBUFFER;
	 nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(lQ,rQ), m));
	 mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), nonZeroBlocks.size()-1));
       }

     }

 }


 /*
 void SpinAdapted::StackWavefunction::FlattenInto (Matrix& C)
 {
   C.ReSize(totalMemory,1);
   DCOPY(totalMemory, data, 1, C.Store(), 1);
 }
 */

 void SpinAdapted::StackWavefunction::CollectFrom (const RowVector& C)
 {
   assert(totalMemory == C.Storage());
   long index = 0;
   for (int lQ = 0; lQ < nrows(); ++lQ)
   for (int rQ = 0; rQ < ncols(); ++rQ)
   if (allowed(lQ, rQ)) {
     int copysize = operator_element(lQ, rQ).Nrows()* operator_element(lQ, rQ).Ncols();
     DCOPY(copysize, C.Store()+(index), 1, operator_element(lQ, rQ).Store(), 1);
     index += copysize;
   }
 }



 void SpinAdapted::StackWavefunction::SaveWavefunctionInfo (const StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num)
 {
   dmrginp.diskwo->start();
   char file [5000];
   int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
   sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
   p1out << "\t\t\t Saving Wavefunction " << file << endl;
   if (mpigetrank() == 0)
     {
       std::ofstream ofs(file, std::ios::binary);
       boost::archive::binary_oarchive save_wave(ofs);
       save_wave << onedot << waveInfo << *waveInfo.leftStateInfo << *(waveInfo.leftStateInfo->leftStateInfo);
       save_wave << *(waveInfo.leftStateInfo->rightStateInfo) << *waveInfo.rightStateInfo;
       if(!onedot)
	 save_wave << *(waveInfo.rightStateInfo->leftStateInfo) << *(waveInfo.rightStateInfo->rightStateInfo);

       this->Save (ofs);
       //save_wave << operatorMatrix;
       ofs.close();
     }
   dmrginp.diskwo->stop();

 }

bool SpinAdapted::StackWavefunction::exists(int state) {
  char file[5000];
  int first = 0, last = 1;
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.load_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", state, ".tmp");
  boost::filesystem::path wavefile(file);
  return boost::filesystem::exists(wavefile);
}

void SpinAdapted::StackWavefunction::CopyState(int from, int to) {
  std::vector<int> sites;
  int new_site, wave_site;
  char file[5000];
  int first = 0, last = 1;
  sprintf (file, "%s", dmrginp.load_prefix().c_str());
  //using namepsace boost::filesystem;
  boost::filesystem::path dir_path(file);

  boost::filesystem::directory_iterator end_itr;
  for ( boost::filesystem::directory_iterator itr( dir_path );
        itr != end_itr;
        ++itr )
    {
      if ( is_directory(itr->status()) )
	  continue;
      else {
	string filename =itr->path().filename().string();
	int len = filename.length();
	string wavesubstrfrom = "."+to_string(from)+".tmp";
	string wavesubstrto = "."+to_string(to)+".tmp";
	std::size_t found = filename.find(wavesubstrfrom);
	if (found != std::string::npos && filename.compare(0,4, "wave")==0) {
	  string outfile = filename;
	  outfile.replace(found, wavesubstrfrom.length(), wavesubstrto.c_str());

	  char filein[5000];
	  sprintf (filein, "%s/%s", dmrginp.load_prefix().c_str() , filename.c_str());
	  boost::filesystem::path frompath(filein);
	  char fileout[5000];
	  sprintf (fileout, "%s/%s", dmrginp.load_prefix().c_str() , outfile.c_str());
	  boost::filesystem::path topath(fileout);

	  copy_file(frompath, topath,boost::filesystem::copy_option::overwrite_if_exists);
	}

	string filename2 = itr->path().filename().string();
	string rotsubstrfrom = ".state"+to_string(from)+".tmp";
	string rotsubstrto = ".state"+to_string(to)+".tmp";
	std::size_t found2 = filename2.find(rotsubstrfrom);
	if (found2 != std::string::npos && filename2.compare(0,8, "Rotation")==0) {
	  string outfile = filename2;
	  outfile.replace(found2, rotsubstrfrom.length(), rotsubstrto.c_str());

	  char filein[5000];
	  sprintf (filein, "%s/%s", dmrginp.load_prefix().c_str() , filename2.c_str());
	  boost::filesystem::path frompath(filein);
	  char fileout[5000];
	  sprintf (fileout, "%s/%s", dmrginp.load_prefix().c_str() , outfile.c_str());
	  boost::filesystem::path topath(fileout);

	  copy_file(frompath, topath,boost::filesystem::copy_option::overwrite_if_exists);
	}

      }
    }


}

void SpinAdapted::StackWavefunction::ChangeLastSite(int newLast, int oldLast, int state) {
  std::vector<int> sites;
  int new_site, wave_site;
  char file[5000];
  int first = 0, last = 1;
  sprintf (file, "%s", dmrginp.load_prefix().c_str());
  //using namepsace boost::filesystem;
  boost::filesystem::path dir_path(file);

  boost::filesystem::directory_iterator end_itr;
  for ( boost::filesystem::directory_iterator itr( dir_path );
        itr != end_itr;
        ++itr )
    {
      if ( is_directory(itr->status()) )
	  continue;
      else {
	string filename =itr->path().filename().string();
	int len = filename.length();
	string wavesubstrfrom = to_string(oldLast)+"."+to_string(mpigetrank())+"."+to_string(state)+".tmp";
	string wavesubstrto   = to_string(newLast)+"."+to_string(mpigetrank())+"."+to_string(state)+".tmp";
	std::size_t found = filename.find(wavesubstrfrom);
	if (found != std::string::npos && filename.compare(0,4, "wave")==0) {
	  string outfile = filename;
	  outfile.replace(found, wavesubstrfrom.length(), wavesubstrto.c_str());

	  char filein[5000];
	  sprintf (filein, "%s/%s", dmrginp.load_prefix().c_str() , filename.c_str());
	  boost::filesystem::path frompath(filein);
	  char fileout[5000];
	  sprintf (fileout, "%s/%s", dmrginp.load_prefix().c_str() , outfile.c_str());
	  boost::filesystem::path topath(fileout);

	  copy_file(frompath, topath,boost::filesystem::copy_option::overwrite_if_exists);
	}

	string filename2 = itr->path().filename().string();
	string rotsubstrfrom = to_string(oldLast)+"."+to_string(mpigetrank())+".state"+to_string(state)+".tmp";
	string rotsubstrto   = to_string(newLast)+"."+to_string(mpigetrank())+".state"+to_string(state)+".tmp";
	std::size_t found2 = filename2.find(rotsubstrfrom);
	if (found2 != std::string::npos && filename2.compare(0,8, "Rotation")==0) {
	  string outfile = filename2;
	  outfile.replace(found2, rotsubstrfrom.length(), rotsubstrto.c_str());

	  char filein[5000];
	  sprintf (filein, "%s/%s", dmrginp.load_prefix().c_str() , filename2.c_str());
	  boost::filesystem::path frompath(filein);
	  char fileout[5000];
	  sprintf (fileout, "%s/%s", dmrginp.load_prefix().c_str() , outfile.c_str());
	  boost::filesystem::path topath(fileout);

	  copy_file(frompath, topath,boost::filesystem::copy_option::overwrite_if_exists);
	}

      }
    }


}

 void SpinAdapted::StackWavefunction::LoadWavefunctionInfo (StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num, bool allocateData)
 {
   dmrginp.diskwi->start();
   char file [5000];
   int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
   sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.load_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
   p1out << "\t\t\t Loading Wavefunction " << file << endl;
   waveInfo.Allocate ();
   if (mpigetrank() == 0)
     {
       std::ifstream ifs(file, std::ios::binary);
       boost::archive::binary_iarchive load_wave(ifs);
       load_wave >> onedot >> waveInfo >> *waveInfo.leftStateInfo >> *(waveInfo.leftStateInfo->leftStateInfo)
		 >> *(waveInfo.leftStateInfo->rightStateInfo) >> *waveInfo.rightStateInfo;
       if(!onedot)
	 load_wave >> *(waveInfo.rightStateInfo->leftStateInfo) >> *(waveInfo.rightStateInfo->rightStateInfo);

       this->Load (ifs, allocateData);
       //load_wave >> operatorMatrix;
       ifs.close();
       allocateWfnOperatorMatrix(); 
     }
   dmrginp.diskwi->stop();
 }

 void SpinAdapted::StackWavefunction::CollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
 {
   //mdebugcheck("before collectquantaalongrows");
   try
     {
       Timer ctimer;

       std::vector<std::pair<std::pair<int, int>, StackMatrix> > nonZeroBlocksbkp = nonZeroBlocks;
       std::map< std::pair<int, int>, int> mapToNonZeroBlocksbkp = mapToNonZeroBlocks;

       StateInfo tmpState = sRow;
       tmpState.CollectQuanta ();

       StackWavefunction tmpOper;
       tmpOper.initialise(deltaQuantum, tmpState, sCol, onedot);

       rowCompressedForm.clear();rowCompressedForm.resize(tmpState.quanta.size(), vector<int>());
       colCompressedForm.clear();colCompressedForm.resize(sCol.quanta.size(), vector<int>());
       nonZeroBlocks.resize(0);
       mapToNonZeroBlocks.clear();


       ObjectMatrix<StackMatrix*> matRef;
       for (int i = 0; i < tmpState.quanta.size (); ++i)
	 for (int j = 0; j < sCol.quanta.size (); ++j) {
	   std::vector<int> dum (1); dum [0] = j;
	   if (tmpOper.allowed(i, j)) {
	     {
	       int rows = tmpState.oldToNewState[i].size ();
	       int cols = dum.size ();
	       matRef.ReSize (rows, cols);
	       for (int x = 0; x < rows; ++x)
		 for (int y = 0; y < cols; ++y)
		   {
		     matRef (x,y) = &nonZeroBlocksbkp[mapToNonZeroBlocksbkp.at(std::pair<int,int>( tmpState.oldToNewState[i][x], dum[y] ))].second;
		     //m (x,y) = &operatorMatrix(oldToNewStateI [i], oldToNewStateJ [j]);
		   }
	     }
	     //OperatorMatrixReference (matRef, tmpState.oldToNewState [i], dum);
	     CatenateProduct (matRef, tmpOper.operator_element(i,j));

	     rowCompressedForm[i].push_back(j);
	     colCompressedForm[j].push_back(i);
	     nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(i,j), tmpOper.operator_element(i,j)));
	     mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(i, j), nonZeroBlocks.size()-1));
	   }
	 }

       allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
       //operatorMatrix = tmpOper.operatorMatrix;

       allocateWfnOperatorMatrix(); 
       copyData(tmpOper);
       //copy(tmpOper.operatorMatrix, operatorMatrix);
       //DCOPY(totalMemory, tmpOper.data, 1, data, 1);

       tmpOper.deallocate();
     }
   catch (Exception)
     {
       Exception::what ();
       abort ();
     }
   //mdebugcheck("after collectquantaalongrows");
 }



 void SpinAdapted::StackWavefunction::UnCollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
 {
   try
     {
       std::vector<std::pair<std::pair<int, int>, StackMatrix> > nonZeroBlocksbkp = nonZeroBlocks;
       std::map< std::pair<int, int>, int> mapToNonZeroBlocksbkp = mapToNonZeroBlocks;

       rowCompressedForm.clear();rowCompressedForm.resize(sRow.unCollectedStateInfo->quanta.size(), vector<int>());
       colCompressedForm.clear();colCompressedForm.resize(sCol.quanta.size(), vector<int>());
       nonZeroBlocks.resize(0);
       mapToNonZeroBlocks.clear();

       StackWavefunction tmpOper;
       tmpOper.initialise(deltaQuantum, *sRow.unCollectedStateInfo, sCol, onedot);
       tmpOper.Clear();

       for (int i = 0; i < sRow.quanta.size (); ++i)
	 {
	   const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	   int firstRow = 0;
	   for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	     {
	       int unCollectedI = oldToNewStateI [iSub];
	       int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	       for (int j = 0; j < sCol.quanta.size (); ++j)
		 if (tmpOper.allowedQuantaMatrix (unCollectedI, j)) {
		   for (int row = 0; row<lastRowSize; row++) 
		     for (int col = 0; col <tmpOper.operator_element(unCollectedI, j).Ncols(); col++) {
		       tmpOper.operator_element(unCollectedI, j)(row+1, col+1) = nonZeroBlocksbkp[mapToNonZeroBlocksbkp.at(std::pair<int,int>(i,j))].second( row + firstRow + 1, col+1);
		       //tmpOper.operator_element(unCollectedI, j)(row+1, col+1) = operatorMatrix(i, j)( row + firstRow + 1, col+1);
		     }


		  rowCompressedForm[unCollectedI].push_back(j);
		  colCompressedForm[j].push_back(unCollectedI);
		  nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(unCollectedI,j), tmpOper.operator_element(unCollectedI, j)));
		  mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(unCollectedI, j), nonZeroBlocks.size()-1));
		}
	      firstRow += lastRowSize;
	    }
	}
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      //operatorMatrix = tmpOper.operatorMatrix;

      allocateWfnOperatorMatrix(); 
      copyData(tmpOper);
      //copy(tmpOper.operatorMatrix, operatorMatrix);
      
      tmpOper.deallocate();

    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
}



void SpinAdapted::StackWavefunction::CollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  //mdebugcheck("before collectquantaalongcolumns");
  try
    {
      std::vector<std::pair<std::pair<int, int>, StackMatrix> > nonZeroBlocksbkp = nonZeroBlocks;
      std::map< std::pair<int, int>, int> mapToNonZeroBlocksbkp = mapToNonZeroBlocks;

      StateInfo tmpState = sCol;
      tmpState.CollectQuanta ();

      StackWavefunction tmpOper;
      tmpOper.initialise (deltaQuantum, sRow, tmpState, onedot);
      tmpOper.Clear();

      rowCompressedForm.clear();rowCompressedForm.resize(sRow.quanta.size(), vector<int>());
      colCompressedForm.clear();colCompressedForm.resize(tmpState.quanta.size(), vector<int>());
      nonZeroBlocks.resize(0);
      mapToNonZeroBlocks.clear();

      ObjectMatrix<StackMatrix*> matRef;
      for (int i = 0; i < sRow.quanta.size (); ++i)
	for (int j = 0; j < tmpState.quanta.size (); ++j)
	  {
	    std::vector<int> dum (1); dum [0] = i;
	    if (tmpOper.allowed(i, j))
	      {
		{
		  int rows = dum.size ();
		  int cols = tmpState.oldToNewState[j].size ();
		  matRef.ReSize (rows, cols);
		  for (int x = 0; x < rows; ++x)
		    for (int y = 0; y < cols; ++y)
		      {
			matRef (x,y) = &nonZeroBlocksbkp[mapToNonZeroBlocksbkp.at(std::pair<int,int>(dum[x],tmpState.oldToNewState[j][y]))].second;
			//m (x,y) = &operatorMatrix(oldToNewStateI [i], oldToNewStateJ [j]);
		      }
		}
		//OperatorMatrixReference (matRef, dum, tmpState.oldToNewState [j]);
		CatenateProduct (matRef, tmpOper.operator_element(i,j));

		rowCompressedForm[i].push_back(j);
		colCompressedForm[j].push_back(i);
		nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(i,j), tmpOper.operator_element(i,j)));
		mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(i, j), nonZeroBlocks.size()-1));
	      }
	  }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      //operatorMatrix = tmpOper.operatorMatrix;

      allocateWfnOperatorMatrix(); 
      copyData(tmpOper);
      //copy(tmpOper.operatorMatrix, operatorMatrix);
      //DCOPY(totalMemory, tmpOper.data, 1, data, 1);

      tmpOper.deallocate();
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongcolumns");
}

void SpinAdapted::StackWavefunction::deepCopy(const StackWavefunction& o)
{
  *this = o;
  data = Stackmem[omprank].allocate(totalMemory);
  DCOPY(totalMemory, o.data, 1, data, 1);
  allocateWfnOperatorMatrix();
}

void SpinAdapted::StackWavefunction::deepClearCopy(const StackWavefunction& o)
{
  *this = o;
  data = Stackmem[omprank].allocate(totalMemory);
  memset(data, 0, totalMemory*sizeof(double));
  allocateWfnOperatorMatrix();
}

void  SpinAdapted::StackWavefunction::UnCollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  std::vector<std::pair<std::pair<int, int>, StackMatrix> > nonZeroBlocksbkp = nonZeroBlocks;
  std::map< std::pair<int, int>, int> mapToNonZeroBlocksbkp = mapToNonZeroBlocks;

  rowCompressedForm.clear();rowCompressedForm.resize(sRow.quanta.size(), vector<int>());
  colCompressedForm.clear();colCompressedForm.resize(sCol.unCollectedStateInfo->quanta.size(), vector<int>());
  nonZeroBlocks.resize(0);
  mapToNonZeroBlocks.clear();

  StackWavefunction tmpOper;
  tmpOper.initialise(deltaQuantum, sRow, *sCol.unCollectedStateInfo, onedot);
  tmpOper.Clear();

  for (int j = 0; j < sRow.quanta.size (); ++j)
  for (int i = 0; i < sCol.quanta.size (); ++i)
  {
    const std::vector<int>& oldToNewStateI = sCol.oldToNewState [i];
    int firstCol = 0;
    for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
    {
      int unCollectedI = oldToNewStateI [iSub];
      int lastColSize = sCol.unCollectedStateInfo->quantaStates [unCollectedI];
      if (tmpOper.allowed(j, unCollectedI)){
	for (int row = 0; row<tmpOper.operator_element(j, unCollectedI).Nrows(); row++) 
	for (int col = 0; col <lastColSize; col++) 
	  tmpOper.operator_element(j, unCollectedI)(row+1, col+1) = nonZeroBlocksbkp[mapToNonZeroBlocksbkp.at(std::pair<int,int>(j,i))].second( row + 1, col+1+firstCol);
	//operatorMatrix(j, i)( row + 1, col+1+firstCol);
	      	      
	rowCompressedForm[j].push_back(unCollectedI);
	colCompressedForm[unCollectedI].push_back(j);
	nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(j, unCollectedI), tmpOper.operator_element(j, unCollectedI)));
	mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(j, unCollectedI), nonZeroBlocks.size()-1));
	
      }
      firstCol += lastColSize;
    }
  }
  allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
  //operatorMatrix = tmpOper.operatorMatrix;
  allocateWfnOperatorMatrix(); 
  copyData(tmpOper);
  //copy(tmpOper.operatorMatrix, operatorMatrix);

  
  tmpOper.deallocate();
}




SpinAdapted::StackWavefunction& SpinAdapted::StackWavefunction::operator+=(const StackWavefunction& other)
{
  for (int i = 0; i < nrows(); ++i)
    for (int j = 0; j < ncols(); ++j)
      if (allowed(i, j))
        {
          assert(other.allowed(i, j));
          MatrixScaleAdd(1., other.operator_element(i, j), operator_element(i, j));
        }
  return *this;
}
