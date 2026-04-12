function [nnbors,nbors] = computecoll(ndim,nlevels,nboxes,laddr,boxsize,centers,iparent,nchild,ichild,iperiod,nnbors,nbors)

% laddr(2,0:nlevels)
% boxsize(0:nlevels)

idx1 = 1;

bs0=boxsize(0+idx1);

mc= 2^ndim;
mnbors=3^ndim;

for i=1:nboxes
  nnbors(i) = 0;
  for j=1:mnbors
    nbors(j,i) = -1;
  end
end

nnbors(1) = 1;
nbors(1,1) = 1;
for ilev = 1:nlevels
  ifirstbox = laddr(1,ilev+idx1);
  ilastbox = laddr(2,ilev+idx1);
  
  for ibox = ifirstbox:ilastbox
    dad = iparent(ibox);
    for i=1:nnbors(dad)
      jbox = nbors(i,dad);
      for j=1:mc
        kbox = ichild(j,jbox);
        if(kbox>0)
          ifnbor=1;
          for k=1:ndim
            dis=abs(centers(k,kbox)-centers(k,ibox));
            if (iperiod==1)
              dp1=bs0-dis;
              if (dp1<dis), dis=dp1; end
            end
            if (dis>1.05*boxsize(ilev+idx1))
              ifnbor=0;
              break
            end
          end
             
          if(ifnbor==1)
            nnbors(ibox) = nnbors(ibox)+1;
            nbors(nnbors(ibox),ibox) = kbox;
          end
        end
      end
    end
  end
end

end