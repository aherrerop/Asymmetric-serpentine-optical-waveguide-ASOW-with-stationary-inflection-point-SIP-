function Min_HD = SIP_Check_Eigenvectors(Teq)
%this fun calculates min HD for 3 vectors out of 6 vectors 
%Teq should be 6x6 matrix

[V,~]=eig(Teq);%calculate eigenvector matrix

C = combnk(1:6,3) ;%get combination of 3 out of 6 

for index=1:length(C)
    
    V1=V(:,C(index,1));
    V2=V(:,C(index,2));
    V3=V(:,C(index,3));
   
    C12=(1-abs(V2'*V1)^2)^0.5;
    C13=(1-abs(V3'*V1)^2)^0.5;
    C23=(1-abs(V3'*V2)^2)^0.5;
    HD(index)=mean(abs([C12 C13 C23]));

end


Min_HD=min(HD);


