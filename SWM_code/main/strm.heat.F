c this code reads in "daily.tmp" from "main.F" and 
c computes streamfunction responses

      include "param.h"  
      include "const1.h" 

      parameter(nv=14)  
      parameter(id0=1)   
      parameter(id1=60-1)
      parameter(id2=1)
      
      complex vor1(mx,jx,nv),div1(mx,jx,nv),
     *        t1(mx,jx,nv),z1(mx,jx,nv),ps1(mx,jx),
     *        vorm1(mx,jx),str1(mx,jx)

      real vorgbar(il,nv),
     * vorg1(ix,il,nv),vorgvar(ix,il,nv) 

      real stream(ix,il,nv),avgstream(ix,il,nv)

      character*75 fname,fname1,fname2  

      call spinit
      call ginit

C...OUTPUT  

      open(1,file=
c     &'strm.tran.tmp',
     &'../output/strm.heat.tmp',
     &  form='unformatted',access='direct',recl=ix*il*kx*4)

C...INPUT 
      open(80,file= 
     &'../temp/outfile1', 
     &  form='unformatted')  

c...skip the first id0 records

      do index=1,id0 
        read(80) vor1,div1,t1,z1,ps1
      end do

c...read in id1*id2 records, calculate id1 times every id2 records

      do 999 index=1,id1     
      print*,'index=',index 

      do k=1,nv  
      do j=1,il
        vorgbar(j,k)=0.0
        do i=1,ix
          vorgvar(i,j,k)=0.0
        end do
      end do
      end do

      read(80) vor1,div1,t1,z1,ps1

      do k=1,nv
        call gridd(vor1(1,1,k),vorg1(1,1,k),1)

        do j=1,il
        do i=1,ix
	   vorgbar(j,k)=vorgbar(j,k)+vorg1(i,j,k)/float(ix)
        end do
        end do
      end do

      do i=1,ix
      do j=1,il
      do k=1,nv
         vorgvar(i,j,k)=vorg1(i,j,k)-vorgbar(j,k)
      end do
      end do
      end do

      do k=1,nv
         call spec(vorgvar(1,1,k),vorm1)
         call invlap(vorm1,str1)
         call gridd(str1,stream(1,1,k),1)
      end do

	write(1,rec=index) stream 

 999  continue      
      close(80)  
   
      stop
      end

