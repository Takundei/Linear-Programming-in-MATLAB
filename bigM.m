function solution=bigM(c,A,b)
[tableau, basis,mR] = makeTableau(c,A,b);
simplexSol = simplexMethod(tableau,basis);

if simplexSol==0
    solution=0;
else 
  x=round(simplexSol,5);
 z=x(end); 
 test=sum(x(end-mR:end-1));
 if test==0
    x=x(1:length(c));
  solution=[x ;z]; 
 else
 solution=1;    %infeasible
 end
end

function [tableau, basis,mR] = makeTableau(c,A,b)
[n,~]=size(A);
%changing to +b since b can be possitive or negative
for i=1:length(b)
    if b(i)<0
        A(i,:)=-1*A(i,:);
        b(i)=-1*b(i);
    end
end
d=ones(n,1);
lA=length(A);
mR=0;
[wA,~]=size(A);
 M=99999999;
 % adding Surplus variables
for i=1:wA
    e=zeros(n,1);
    e(i)=1;
    count=0;
    for j=1:lA
    if A(i,j)~=1 || sum(abs(A(:,j)))~=1
        count=count +1;
    end
    end
%Adding the Big M and     
    if count==lA
        mR=mR+1;
        A(:,lA+mR)=e;
        bigM(mR)=M;
    end
end
E=[A,b];
p=lA-length(c);%fills the space from where c ends to the first M with zeros
r1=[-c',zeros(1,p),bigM,0];
tableau=[r1;E];
T=tableau(2:end,:);
[m,l]=size(T);
q=1;
% changing to canonical form, removing the M's above the basic var columns
    for i=1:mR
    C=tableau(:,lA+i);
    [R,~]=find(C==1);
    tableau(1,:)=tableau(1,:)-tableau(1,lA+i)*tableau(R,:);
    end
%forming the basis    
for i=1:l
if T(:,i)\d==1
    basis(q)=i;
    q=q+1;
else
    continue
end
end
basis=[basis]';
end
function simplexSol = simplexMethod(tableau,basis,c)
pivotElement = getPivot(tableau);
if pivotElement==0 
    simplexSol='Infeasible LP';
else
    while pivotElement~=0
        if pivotElement==1
            simplexSol=0;
             break
        else 
            [tableau, basis] = updateTableau(tableau, basis, pivotElement);
        end
        pivotElement = getPivot(tableau);
 [x,z] = readSolution(tableau,basis);
  simplexSol=[x ;z]; 
 end

 %For the multiple solutions...
%  [s,~]=size(solution);
%  if s>1
%   Alt=0;
%  of=tableau(1,:);
%  for i=1:length(basis)
%      of(basis(i))=1;
%  end
%  for i=1:length(z)
%      if of(i)==0
%          Alt=i;
%      end
%  end
%  if Alt>0
%      pivotRow=getPivotRow(tableau, Alt);
%      pivotElement = [pivotRow+1;Alt];
%      [tableau, basis] = updateTableau(tableau, basis, pivotElement);
%  end
%  [tu,u] = readSolution(tableau,basis);
%  solution(:,2)=[tu;u];
%  end
end
%  function pivotRow= getPivotRow(tableau, pivotCol)
% pivotRow = 0;
% minFound = Inf;
% b=tableau(:,end);
% a=tableau(:,pivotCol);
% for i=2:length(b)
%     if b(i)>0 && a(i)/b(i)<minFound && a(i)/b(i)>0
%         minFound = a(i)/b(i);
%         pivotRow = i-1;
%     end
% end
% end
             
function pivotElement = getPivot(tableau)
pivotCol=getPivotColumn(tableau);

if getPivotColumn(tableau)==0
    pivotElement=0;
elseif getPivotRow(tableau, pivotCol)==0
    pivotElement=1;
else
    pivotRow=getPivotRow(tableau, pivotCol);
    pivotElement=[pivotRow;pivotCol];
end
%Getting our pivot Column
function pivotCol = getPivotColumn(tableau)
t=tableau(1,1:end-1);
min=0;
pivotCol=0;
for i=1:length(t)
    if t(i)< 0 && t(i)<min
         min=t(i);
            pivotCol=i;
    end
end
end

%Getting our pivot row
function pivotRow= getPivotRow(tableau, pivotCol)
pivotRow = 0;
minFound = 999999;
b=tableau(:,end);
a=tableau(:,pivotCol);
for i=2:length(b)
    if b(i)>0 && b(i)/a(i)>0 && b(i)/a(i)<=minFound
        minFound = b(i)/a(i);
        pivotRow = i;
    end
end
end

end

function [newSimplex, newBasis] = updateTableau(tableau, basis, pivotElement)

    a=pivotElement(1); %+1
    b=pivotElement(2);
    tableau(a,:)=tableau(a,:)/tableau(a,b);
    tableauPivotRow=tableau(a,:);
    [n,~] = size(tableau);
    
newBasis = updateBasis(basis, pivotElement);    

    for i=1:n
        if i~=a
            tableau(i,:)=updateTableauRow(tableauPivotRow,tableau(i,:), pivotElement(2));
        end
    end
newSimplex=tableau;
    function newBasis = updateBasis(basis, pivotElem)
        r=pivotElem(1);
        basis(r)=pivotElem(2);
        newBasis=basis;
    end

    function newRow = updateTableauRow(tableauPivotRow,tableauRow, pivotCol)

        tableauPivotRow=tableauPivotRow/tableauPivotRow(pivotCol);
        pivot=tableauRow(pivotCol)/tableauPivotRow(pivotCol);
        newRow=tableauRow-pivot*tableauPivotRow;

    end
end

function [x,z] = readSolution(tableau,basis)
x=zeros(1,length(tableau)-1)';
for i=1:length(basis)
    n=basis(i);
x(n)=tableau(:,n)\tableau(:,end);
end
z=tableau(1,end);
end
 end
end