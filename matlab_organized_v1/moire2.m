clear all;
clc


a=1;W=1;Lx=10;
t=1;delta0=a/log(10);
al=0;by=t;Ub=0;db=0;dt=0;




mb=5;nb=5;
mt=5;nt=4;
%%%%%%%%%%%%%%%%%%%%%%%
tau1=mb/nb;tau2=mt/nt;
a_strain1=a*tau1;
a_strain2=a*tau2;
t1=t*exp(-(a_strain1-a)/delta0);
t2=t*exp(-(a_strain2-a)/delta0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ab1=[a,0];
ab2=[0,a];
at1=a*[nb/nt,0];
at2=[0,a];
[Rbx,Rby] = B_Lattice_UnitCell(ab2,W);

[Rtx,Rty] = T_Lattice_UnitCell(at2,W);




[RBx,RBy] = B_Lattice_SuperCell(ab1,nb,Rbx,Rby);
[RBX0,RBY0] = B_Lattice_AllCell(ab1,Lx,Rbx,Rby);
RBX=RBX0-RBX0(ceil(Lx/2));RBY=RBY0;
[RTx,RTy] = T_Lattice_SuperCell(at1,nt,Rtx,Rty);
[RTX0,RTY0] = T_Lattice_AllCell(at1,Lx,Rtx,Rty);
RTX=RTX0-RTX0(ceil(Lx/2));RTY=RTY0+at2(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Hb,alphab,betab] = Hamiltonian_Square(al,t1,by,Ub,W,Lx);
[Ht,alphat,betat] = Hamiltonian_Square(al,t2,by,Ub,W,Lx);
[H_inter] = Ugg_Inter(t,RTX,RBX,RTY,RBY);

H=[Hb,H_inter';
    H_inter,Ht];
[ES,EV]=eig(H);
DEV=diag(EV);
[DEV2,index]=sort(DEV,'ComparisonMethod','real');%%energy
REV=real(DEV2);
IEV=imag(DEV2);
ES=ES(:,index);




figure;
scatter(RBX(:),RBY(:)...
        ,100,'filled','r',...
'MarkerFaceAlpha',0.8,...
'MarkerEdgeAlpha',0.8);

hold on;
scatter(RTX(:),RTY(:)...
        ,50,'filled','b',...
'MarkerFaceAlpha',0.8,...
'MarkerEdgeAlpha',0.8);
% 
%     hold off;
% set(gca,'xaxislocation','origin','yaxislocation','origin') % set origin position
% %     xlim([-8,8]);ylim([-8,8]);
% 
% % title(['W_1=',num2str(W1),', L_1=',num2str(NL1),...
% %      ', W_2=',num2str(W2),...
% %      ', \tau_1=',num2str(tau1),', \tau_2=',num2str(tau2)]);
hold on
% quiver(0,0,ab1(1),ab1(2),1,'M','linewidth',3,'MaxHeadSize',1.2,'AutoScaleFactor',1,'AutoScale','off');
% hold on
% quiver(0,0,ab2(1),ab2(2),1,'M','linewidth',3,'MaxHeadSize',1.2,'AutoScaleFactor',1,'AutoScale','off');
% hold on
% quiver(0,0,at1(1),at1(2),1,'g','linewidth',3,'MaxHeadSize',1.2,'AutoScaleFactor',1,'AutoScale','off');
% hold on
% quiver(0,0,at2(1),at2(2),1,'g','linewidth',3,'MaxHeadSize',1.2,'AutoScaleFactor',1,'AutoScale','off');
hold off;
axis equal;
axis off
grid on;
set(gca,'FontSize',15,'Fontname', 'Times New Roman','linewidth',2);
ax = gca;
ax.LineWidth = 2;
ax.XColor = 'k';
ax.YColor = 'k';
set (gcf,'Position',[100,100,400,400]);





function [RBx0,RBy0] = B_Lattice_UnitCell(ab2,W)

    for wi=1:W
    RBy0(wi)=(wi-1)*ab2(2);
    RBx0(wi)=(wi-1)*ab2(1);
    end
    

end



function [RTx0,RTy0] = T_Lattice_UnitCell(at2,W)

    for wi=1:W
    RTy0(wi)=(wi-1)*at2(2);
    RTx0(wi)=(wi-1)*at2(1);
    end
    

end

function [RBx,RBy] = B_Lattice_SuperCell(ab1,nb,Rbx,Rby)

for SCi=1:nb
    RBx(:,SCi)=Rbx+(SCi-1)*ab1(1);
    RBy(:,SCi)=Rby+(SCi-1)*ab1(2);
end


end


function [RBX,RBY] = B_Lattice_AllCell(ab1,Lx,Rbx,Rby)


    NB=[0:1:Lx-1];
    for nbi0=1:1*Lx
        nbi=NB(nbi0);
        
        RBX(:,nbi0)=Rbx+(nbi)*ab1(1);
        RBY(:,nbi0)=Rby+(nbi)*ab1(2);

        
    end




end


function [RTx,RTy] = T_Lattice_SuperCell(at1,nt,Rtx,Rty)

for SCi=1:nt
    RTx(:,SCi)=Rtx+(SCi-1)*at1(1);
    RTy(:,SCi)=Rty+(SCi-1)*at1(2);
end


end


function [RTX,RTY] = T_Lattice_AllCell(at1,Lx,Rtx,Rty)


    NT=[0:1:Lx-1];
    for nti0=1:1*Lx
        nti=NT(nti0);
        
        RTX(:,nti0)=Rtx+(nti)*at1(1);
        RTY(:,nti0)=Rty+(nti)*at1(2);

        
    end




end




function [H,alpha,beta] = Hamiltonian_Square(al,bx,by,Us,W,NL)


alpha=kron(diag(ones(1,W)),al+Us)+...
    kron(diag(ones(1,W-1),1),by')+...
    kron(diag(ones(1,W-1),-1),by);
beta=(bx)*diag(ones(1,W));

H=kron(diag(ones(1,NL)),alpha)+...
    kron(diag(ones(1,NL-1),1),beta')+...
    kron(diag(ones(1,NL-1),-1),beta);
end



function [Ugg0] = Ugg_Inter(t,Rtx,Rbx,Rty,Rby)
V0pi=t;
% V0pi=0;V0si=0;
a=1;delta0=a/log(10);

    for Rti=1:length(Rtx)
        for Rbi1=1:length(Rbx)
            
                r=sqrt((Rtx(Rti)-Rbx(Rbi1))^2+(Rty(Rti)-Rby(Rbi1))^2);
                Ugg0(Rti,Rbi1)=(V0pi*exp(-(r-a)/delta0));
            
            
    
        end
        
        
    end

end
