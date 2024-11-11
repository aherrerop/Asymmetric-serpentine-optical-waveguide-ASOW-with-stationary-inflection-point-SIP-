function [ee,e1,e2] = SIP_Check_Vectors2(Teq)

[V,~]=eig(Teq);

C = combnk(1:6,3);
for i=1:length(C)
    
    I=[C(i,1) C(i,2) C(i,3)];
   
    V1=V(:,I(1));
    V2=V(:,I(2));
    V3=V(:,I(3));
    C12=abs(acos(abs(sum(V1.*conj(V2)))));
    C23=abs(acos(abs(sum(V2.*conj(V3)))));
    C31=abs(acos(abs(sum(V3.*conj(V1)))));
    e1=(C12^2+C23^2+C31^2)^0.5;
  
    VV2=V;
    VV2(:,I)=[];
    V1=VV2(:,1);
    V2=VV2(:,2);
    V3=VV2(:,3);
    C12=acos(abs(sum(V1.*conj(V2))));
    C23=acos(abs(sum(V2.*conj(V3))));
    C31=acos(abs(sum(V3.*conj(V1))));
    e2=(C12^2+C23^2+C31^2)^0.5;
    
     e(i)=(e1^2+e2^2)^0.5;
    
end

ee=min(e);
