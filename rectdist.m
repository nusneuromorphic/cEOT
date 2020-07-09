function norm_dist=rectdist(R1,R2)
%normalized distance between region R1 and R2
overlap=rectint(R1,R2);
if overlap>0
    norm_dist=0;
else
    centre1=[R1(1)+R1(3)/2 R1(2)+R1(4)/2];
    centre2=[R2(1)+R2(3)/2 R2(2)+R2(4)/2];
    distance=centre1-centre2;
    norm_dist=abs(distance(1))/(R1(3)+R2(3))+abs(distance(2))/(R1(4)+R2(4));
end