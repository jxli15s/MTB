clear all;
clc


a=1;W=1;
mb=5;nb=5;
mt=5;nt=4;

ab1=[a,0];
ab2=[0,a];
at1=a*[nb/nt,0];
at2=[0,a];
[Rbx,Rby] = B_Lattice_UnitCell(ab2,W);

[Rtx,Rty] = T_Lattice_UnitCell(at2,W);




[RBx,RBy] = B_Lattice_SuperCell(ab1,nb,Rbx,Rby);
[RBX,RBY] = B_Lattice_AllCell(ab1,nb,Rbx,Rby);

[RTx,RTy] = T_Lattice_SuperCell(at1,nt,Rtx,Rty);
[RTX,RTY] = T_Lattice_AllCell(at1,nt,Rtx,Rty);
RTY=RTY+a;
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


function [RBX,RBY] = B_Lattice_AllCell(ab1,nb,Rbx,Rby)


    NB=[-nb:1:-1, 0:1:nb-1,nb:1:2*nb-1];
    for nbi0=1:3*nb
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


function [RTX,RTY] = T_Lattice_AllCell(at1,nt,Rtx,Rty)


    NT=[-nt:1:-1, 0:1:nt-1,nt:1:2*nt-1];
    for nti0=1:3*nt
        nti=NT(nti0);
        
        RTX(:,nti0)=Rtx+(nti)*at1(1);
        RTY(:,nti0)=Rty+(nti)*at1(2);

        
    end




end
