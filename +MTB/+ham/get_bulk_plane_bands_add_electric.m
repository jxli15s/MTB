function Energy=get_bulk_plane_bands_add_electric(hamiltonian,hopping_r,wpos,Electric_field_in_evpA,nbands,nrpts,nk,a,b)
    kTrans=@(x)x(1)*b(1,:)+x(2)*b(2,:)+x(3)*b(3,:); %This can be done b'*k'
    x=linspace(-0.5,0.5,nk);
    y=linspace(-0.5,0.5,nk);
    [kx,ky]=meshgrid(x,y);
    kx=reshape(kx,[],1);
    ky=reshape(ky,[],1);
    kpoints=[kx,ky];
    kpoints(:,3)=0;
    kpoints=kpoints*b;

    %%solve the Hamiltonian(k)
    Ek=zeros(nbands,length(kpoints));
    parfor i=1:length(kpoints)
        %fprintf("Processing on %d of %d KPOINTS\n", i,length(kpoints));
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a,wpos,Electric_field_in_evpA);
        [V,D]=eig(hk)
        [E,ind]=sort(diag(D))
        Psik=V(:,ind)
        sort(real(eig(hk)));
        Energy(:,i)=vals;
    end
end

function 

function hk=get_hk(obj,kpoint)
    hk=zeros(size(obj.ham,1),size(obj.ham,1));
    nrpts=size(obj.ham,3);
    for j=1:nrpts
        temp=obj.ham(:,:,j)*...
        exp(1j*dot(kpoint,(obj.hopr(j,:)*obj.a))); 
        hk=hk+temp;
    end
end