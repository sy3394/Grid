/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/SharedMemory.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/GridCore.h>

#ifdef GRID_NVCC
#include <cuda_runtime_api.h>
#endif

NAMESPACE_BEGIN(Grid); 

/*Construct from an MPI communicator*/
void GlobalSharedMemory::Init(Grid_MPI_Comm comm)
{
  assert(_ShmSetup==0);
  WorldComm = comm;
  MPI_Comm_rank(WorldComm,&WorldRank);
  MPI_Comm_size(WorldComm,&WorldSize);
  // WorldComm, WorldSize, WorldRank

  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory
  /////////////////////////////////////////////////////////////////////
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&WorldShmComm);
  MPI_Comm_rank(WorldShmComm     ,&WorldShmRank);
  MPI_Comm_size(WorldShmComm     ,&WorldShmSize);

  std::cout << " GRID initialisation: World communicator of size " <<WorldSize << std::endl;  
  std::cout << " GRID initialisation: Node  communicator of size " <<WorldShmSize << std::endl;
  // WorldShmComm, WorldShmSize, WorldShmRank

  // WorldNodes
  WorldNodes = WorldSize/WorldShmSize;
  assert( (WorldNodes * WorldShmSize) == WorldSize );

  // FIXME: Check all WorldShmSize are the same ?

  /////////////////////////////////////////////////////////////////////
  // find world ranks in our SHM group (i.e. which ranks are on our node)
  /////////////////////////////////////////////////////////////////////
  MPI_Group WorldGroup, ShmGroup;
  MPI_Comm_group (WorldComm, &WorldGroup); 
  MPI_Comm_group (WorldShmComm, &ShmGroup);

  std::vector<int> world_ranks(WorldSize);   for(int r=0;r<WorldSize;r++) world_ranks[r]=r;

  WorldShmRanks.resize(WorldSize); 
  MPI_Group_translate_ranks (WorldGroup,WorldSize,&world_ranks[0],ShmGroup, &WorldShmRanks[0]); 

  ///////////////////////////////////////////////////////////////////
  // Identify who is in my group and nominate the leader
  ///////////////////////////////////////////////////////////////////
  int g=0;
  std::vector<int> MyGroup;
  MyGroup.resize(WorldShmSize);
  for(int rank=0;rank<WorldSize;rank++){
    if(WorldShmRanks[rank]!=MPI_UNDEFINED){
      assert(g<WorldShmSize);
      MyGroup[g++] = rank;
    }
  }
  
  std::sort(MyGroup.begin(),MyGroup.end(),std::less<int>());
  int myleader = MyGroup[0];
  
  std::vector<int> leaders_1hot(WorldSize,0);
  std::vector<int> leaders_group(WorldNodes,0);
  leaders_1hot [ myleader ] = 1;
    
  ///////////////////////////////////////////////////////////////////
  // global sum leaders over comm world
  ///////////////////////////////////////////////////////////////////
  int ierr=MPI_Allreduce(MPI_IN_PLACE,&leaders_1hot[0],WorldSize,MPI_INT,MPI_SUM,WorldComm);
  assert(ierr==0);

  ///////////////////////////////////////////////////////////////////
  // find the group leaders world rank
  ///////////////////////////////////////////////////////////////////
  int group=0;
  for(int l=0;l<WorldSize;l++){
    if(leaders_1hot[l]){
      leaders_group[group++] = l;
    }
  }

  ///////////////////////////////////////////////////////////////////
  // Identify the node of the group in which I (and my leader) live
  ///////////////////////////////////////////////////////////////////
  WorldNode=-1;
  for(int g=0;g<WorldNodes;g++){
    if (myleader == leaders_group[g]){
      WorldNode=g;
    }
  }
  assert(WorldNode!=-1);
  _ShmSetup=1;
}

void GlobalSharedMemory::OptimalCommunicator(const Coordinate &processors,Grid_MPI_Comm & optimal_comm)
{
  ////////////////////////////////////////////////////////////////
  // Assert power of two shm_size.
  ////////////////////////////////////////////////////////////////
  int log2size = -1;
  for(int i=0;i<=MAXLOG2RANKSPERNODE;i++){  
    if ( (0x1<<i) == WorldShmSize ) {
      log2size = i;
      break;
    }
  }
  assert(log2size != -1);

  ////////////////////////////////////////////////////////////////
  // Identify subblock of ranks on node spreading across dims
  // in a maximally symmetrical way
  ////////////////////////////////////////////////////////////////
  int ndimension              = processors.size();
  Coordinate processor_coor(ndimension);
  Coordinate WorldDims = processors; Coordinate ShmDims(ndimension,1);  Coordinate NodeDims (ndimension);
  Coordinate ShmCoor(ndimension);    Coordinate NodeCoor(ndimension);   Coordinate WorldCoor(ndimension);
  int dim = 0;
  for(int l2=0;l2<log2size;l2++){
    while ( (WorldDims[dim] / ShmDims[dim]) <= 1 ) dim=(dim+1)%ndimension;
    ShmDims[dim]*=2;
    dim=(dim+1)%ndimension;
  }

  ////////////////////////////////////////////////////////////////
  // Establish torus of processes and nodes with sub-blockings
  ////////////////////////////////////////////////////////////////
  for(int d=0;d<ndimension;d++){
    NodeDims[d] = WorldDims[d]/ShmDims[d];
  }

  ////////////////////////////////////////////////////////////////
  // Check processor counts match
  ////////////////////////////////////////////////////////////////
  int Nprocessors=1;
  for(int i=0;i<ndimension;i++){
    Nprocessors*=processors[i];
  }
  assert(WorldSize==Nprocessors);

  ////////////////////////////////////////////////////////////////
  // Establish mapping between lexico physics coord and WorldRank
  ////////////////////////////////////////////////////////////////
  int rank;

  Lexicographic::CoorFromIndexReversed(NodeCoor,WorldNode   ,NodeDims);
  Lexicographic::CoorFromIndexReversed(ShmCoor ,WorldShmRank,ShmDims);
  for(int d=0;d<ndimension;d++) WorldCoor[d] = NodeCoor[d]*ShmDims[d]+ShmCoor[d];
  Lexicographic::IndexFromCoorReversed(WorldCoor,rank,WorldDims);

  /////////////////////////////////////////////////////////////////
  // Build the new communicator
  /////////////////////////////////////////////////////////////////
  int ierr= MPI_Comm_split(WorldComm,0,rank,&optimal_comm);
  assert(ierr==0);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Hugetlbfs mapping intended
////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GRID_NVCC
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  void * ShmCommBuf ; 
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the pointer array for shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO/FIXME : NOT ALL NVLINK BOARDS have full Peer to peer connectivity.
  // The annoyance is that they have partial peer 2 peer. This occurs on the 8 GPU blades.
  // e.g. DGX1, supermicro board, 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  cudaDeviceGetP2PAttribute(&perfRank, cudaDevP2PAttrPerformanceRank, device1, device2);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Each MPI rank should allocate our own buffer
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto err =  cudaMalloc(&ShmCommBuf, bytes);
  if ( err !=  cudaSuccess) {
    std::cerr << " SharedMemoryMPI.cc cudaMallocManaged failed for " << bytes<<" bytes " <<cudaGetErrorString(err)<< std::endl;
    exit(EXIT_FAILURE);  
  }
  if (ShmCommBuf == (void *)NULL ) {
    std::cerr << " SharedMemoryMPI.cc cudaMallocManaged failed NULL pointer for " << bytes<<" bytes " << std::endl;
    exit(EXIT_FAILURE);  
  }
  SharedMemoryZero(ShmCommBuf,bytes);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over ranks/gpu's on our node
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int r=0;r<WorldShmSize;r++){
    
    //////////////////////////////////////////////////
    // If it is me, pass around the IPC access key
    //////////////////////////////////////////////////
    cudaIpcMemHandle_t handle;

    if ( r==WorldShmRank ) { 
      err = cudaIpcGetMemHandle(&handle,ShmCommBuf);
      if ( err !=  cudaSuccess) {
	std::cerr << " SharedMemoryMPI.cc cudaIpcGetMemHandle failed for rank" << r <<" "<<cudaGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
    //////////////////////////////////////////////////
    // Share this IPC handle across the Shm Comm
    //////////////////////////////////////////////////
    { 
      int ierr=MPI_Bcast(&handle,
			 sizeof(handle),
			 MPI_BYTE,
			 r,
			 WorldShmComm);
      assert(ierr==0);
    }
    
    ///////////////////////////////////////////////////////////////
    // If I am not the source, overwrite thisBuf with remote buffer
    ///////////////////////////////////////////////////////////////
    void * thisBuf = ShmCommBuf;
    if ( r!=WorldShmRank ) { 
      err = cudaIpcOpenMemHandle(&thisBuf,handle,cudaIpcMemLazyEnablePeerAccess);
      if ( err !=  cudaSuccess) {
	std::cerr << " SharedMemoryMPI.cc cudaIpcOpenMemHandle failed for rank" << r <<" "<<cudaGetErrorString(err)<< std::endl;
	exit(EXIT_FAILURE);
      }
    }
    ///////////////////////////////////////////////////////////////
    // Save a copy of the device buffers
    ///////////////////////////////////////////////////////////////
    WorldShmCommBufs[r] = thisBuf;
  }

  _ShmAllocBytes=bytes;
  _ShmAlloc=1;
}
#else 
#ifdef GRID_MPI3_SHMMMAP
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // allocate the shared windows for our group
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Hugetlbf and others map filesystems as mappable huge pages
  ////////////////////////////////////////////////////////////////////////////////////////////
  char shm_name [NAME_MAX];
  for(int r=0;r<WorldShmSize;r++){
    
    sprintf(shm_name,GRID_SHM_PATH "/Grid_mpi3_shm_%d_%d",WorldNode,r);
    int fd=open(shm_name,O_RDWR|O_CREAT,0666);
    if ( fd == -1) { 
      printf("open %s failed\n",shm_name);
      perror("open hugetlbfs");
      exit(0);
    }
    int mmap_flag = MAP_SHARED ;
#ifdef MAP_POPULATE    
    mmap_flag|=MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
    if ( flags ) mmap_flag |= MAP_HUGETLB;
#endif
    void *ptr = (void *) mmap(NULL, bytes, PROT_READ | PROT_WRITE, mmap_flag,fd, 0); 
    if ( ptr == (void *)MAP_FAILED ) {    
      printf("mmap %s failed\n",shm_name);
      perror("failed mmap");      assert(0);    
    }
    assert(((uint64_t)ptr&0x3F)==0);
    close(fd);
    WorldShmCommBufs[r] =ptr;
  }
  _ShmAlloc=1;
  _ShmAllocBytes  = bytes;
};
#endif // MMAP

#ifdef GRID_MPI3_SHMOPEN
////////////////////////////////////////////////////////////////////////////////////////////
// POSIX SHMOPEN ; as far as I know Linux does not allow EXPLICIT HugePages with this case
// tmpfs (Larry Meadows says) does not support explicit huge page, and this is used for 
// the posix shm virtual file system
////////////////////////////////////////////////////////////////////////////////////////////
void GlobalSharedMemory::SharedMemoryAllocate(uint64_t bytes, int flags)
{ 
  assert(_ShmSetup==1);
  assert(_ShmAlloc==0); 
  MPI_Barrier(WorldShmComm);
  WorldShmCommBufs.resize(WorldShmSize);

  char shm_name [NAME_MAX];
  if ( WorldShmRank == 0 ) {
    for(int r=0;r<WorldShmSize;r++){
	
      size_t size = bytes;
      
      sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",WorldNode,r);
      
      shm_unlink(shm_name);
      int fd=shm_open(shm_name,O_RDWR|O_CREAT,0666);
      if ( fd < 0 ) {	perror("failed shm_open");	assert(0);      }
      ftruncate(fd, size);
	
      int mmap_flag = MAP_SHARED;
#ifdef MAP_POPULATE 
      mmap_flag |= MAP_POPULATE;
#endif
#ifdef MAP_HUGETLB
      if (flags) mmap_flag |= MAP_HUGETLB;
#endif
      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, mmap_flag, fd, 0);
      
      if ( ptr == (void * )MAP_FAILED ) {       perror("failed mmap");      assert(0);    }
      assert(((uint64_t)ptr&0x3F)==0);
      
      WorldShmCommBufs[r] =ptr;
      close(fd);
    }
  }

  MPI_Barrier(WorldShmComm);
  
  if ( WorldShmRank != 0 ) { 
    for(int r=0;r<WorldShmSize;r++){

      size_t size = bytes ;
      
      sprintf(shm_name,"/Grid_mpi3_shm_%d_%d",WorldNode,r);
      
      int fd=shm_open(shm_name,O_RDWR,0666);
      if ( fd<0 ) {	perror("failed shm_open");	assert(0);      }
      
      void * ptr =  mmap(NULL,size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if ( ptr == MAP_FAILED ) {       perror("failed mmap");      assert(0);    }
      assert(((uint64_t)ptr&0x3F)==0);
      WorldShmCommBufs[r] =ptr;

      close(fd);
    }
  }
  _ShmAlloc=1;
  _ShmAllocBytes = bytes;
}
#endif
#endif // End NVCC case for GPU device buffers

/////////////////////////////////////////////////////////////////////////
// Routines accessing shared memory should route through for GPU safety
/////////////////////////////////////////////////////////////////////////
void GlobalSharedMemory::SharedMemoryZero(void *dest,size_t bytes)
{
#ifdef GRID_NVCC
  cudaMemset(dest,0,bytes);
#else
  bzero(dest,bytes);
#endif
}
void GlobalSharedMemory::SharedMemoryCopy(void *dest,const void *src,size_t bytes)
{
#ifdef GRID_NVCC
  cudaMemcpy(dest,src,bytes,cudaMemcpyDefault);
#else   
  bcopy(src,dest,bytes);
#endif
}
////////////////////////////////////////////////////////
// Global shared functionality finished
// Now move to per communicator functionality
////////////////////////////////////////////////////////
void SharedMemory::SetCommunicator(Grid_MPI_Comm comm)
{
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  ShmRanks.resize(size);

  /////////////////////////////////////////////////////////////////////
  // Split into groups that can share memory
  /////////////////////////////////////////////////////////////////////
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,&ShmComm);
  MPI_Comm_rank(ShmComm     ,&ShmRank);
  MPI_Comm_size(ShmComm     ,&ShmSize);
  ShmCommBufs.resize(ShmSize);

  //////////////////////////////////////////////////////////////////////
  // Map ShmRank to WorldShmRank and use the right buffer
  //////////////////////////////////////////////////////////////////////
  assert (GlobalSharedMemory::ShmAlloc()==1);
  heap_size = GlobalSharedMemory::ShmAllocBytes();
  for(int r=0;r<ShmSize;r++){

    uint32_t wsr = (r==ShmRank) ? GlobalSharedMemory::WorldShmRank : 0 ;

    MPI_Allreduce(MPI_IN_PLACE,&wsr,1,MPI_UINT32_T,MPI_SUM,comm);

    ShmCommBufs[r] = GlobalSharedMemory::WorldShmCommBufs[wsr];
  }
  ShmBufferFreeAll();

  /////////////////////////////////////////////////////////////////////
  // find comm ranks in our SHM group (i.e. which ranks are on our node)
  /////////////////////////////////////////////////////////////////////
  MPI_Group FullGroup, ShmGroup;
  MPI_Comm_group (comm   , &FullGroup); 
  MPI_Comm_group (ShmComm, &ShmGroup);

  std::vector<int> ranks(size);   for(int r=0;r<size;r++) ranks[r]=r;
  MPI_Group_translate_ranks (FullGroup,size,&ranks[0],ShmGroup, &ShmRanks[0]); 

  SharedMemoryTest();
}
//////////////////////////////////////////////////////////////////
// On node barrier
//////////////////////////////////////////////////////////////////
void SharedMemory::ShmBarrier(void)
{
  MPI_Barrier  (ShmComm);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test the shared memory is working
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void SharedMemory::SharedMemoryTest(void)
{
  ShmBarrier();
  uint64_t check[3];
  uint64_t magic = 0x5A5A5A;
  if ( ShmRank == 0 ) {
    for(uint64_t r=0;r<ShmSize;r++){
       check[0]=GlobalSharedMemory::WorldNode;
       check[1]=r;
       check[2]=magic;
       GlobalSharedMemory::SharedMemoryCopy( ShmCommBufs[r], check, 3*sizeof(uint64_t));
    }
  }
  ShmBarrier();
  for(uint64_t r=0;r<ShmSize;r++){
    ShmBarrier();
    GlobalSharedMemory::SharedMemoryCopy(check,ShmCommBufs[r], 3*sizeof(uint64_t));
    ShmBarrier();
    assert(check[0]==GlobalSharedMemory::WorldNode);
    assert(check[1]==r);
    assert(check[2]==magic);
    ShmBarrier();
  }
}

void *SharedMemory::ShmBuffer(int rank)
{
  int gpeer = ShmRanks[rank];
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    return ShmCommBufs[gpeer];
  }
}
void *SharedMemory::ShmBufferTranslate(int rank,void * local_p)
{
  int gpeer = ShmRanks[rank];
  assert(gpeer!=ShmRank); // never send to self
  if (gpeer == MPI_UNDEFINED){
    return NULL;
  } else { 
    uint64_t offset = (uint64_t)local_p - (uint64_t)ShmCommBufs[ShmRank];
    uint64_t remote = (uint64_t)ShmCommBufs[gpeer]+offset;
    return (void *) remote;
  }
}

NAMESPACE_END(Grid); 

