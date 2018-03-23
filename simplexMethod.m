 function solution = simplexMethod(tableau,basis)

pivotElement = getPivot(tableau);
if pivotElement==0
    solution='Infeasible LP';
else
    while pivotElement~=0
        if pivotElement==1
            solution=0;
             break
        else 
            [tableau, basis] = updateTableau(tableau, basis, pivotElement)
        end
        pivotElement = getPivot(tableau);
 [x,z] = readSolution(tableau,basis);
 solution=[x ;z];
    end

 %For the multiple solutions...
 [s,~]=size(solution);
 if s>1
  Alt=0;
 of=tableau(1,:);
 for i=1:length(basis)
     of(basis(i))=1;
 end
 for i=1:length(z)
     if of(i)==0
         Alt=i;
     end
 end
 if Alt>0
     pivotRow=getPivotRow(tableau, Alt);
     pivotElement = [pivotRow+1;Alt];
     [tableau, basis] = updateTableau(tableau, basis, pivotElement);
 end
 [tu,u] = readSolution(tableau,basis);
 solution(:,2)=[tu;u];
 end
end
 function pivotRow= getPivotRow(tableau, pivotCol)
pivotRow = 0;
minFound = Inf;
b=tableau(:,end);
a=tableau(:,pivotCol);
for i=2:length(b)
    if b(i)>0 && a(i)/b(i)<minFound && a(i)/b(i)>0
        minFound = a(i)/b(i);
        pivotRow = i-1;
    end
end
end
             
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
minFound = Inf;
b=tableau(:,end);
a=tableau(:,pivotCol);
for i=2:length(b)
    if b(i)>0 && a(i)/b(i)<minFound && a(i)/b(i)>0
        minFound = a(i)/b(i);
        pivotRow = i;
    end
end
end
end

function [newSimplex, newBasis] = updateTableau(tableau, basis, pivotElement)

    a=pivotElement(1)+1;
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
