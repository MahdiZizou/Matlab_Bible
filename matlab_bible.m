%cd('copy paste directory') to change directory

%type ver in console, to see what packages are available for you. 

%script is known as mfile: you can right click on mfile in current folder
%window and run mfiles
%you can type your mfile name in script window and press enter to run it!

%Ctrl+R for multi selected lines commenting
%Ctrl+T for multi selected lines uncommenting
 

edit %this command creat a new script

%shortcuts: ctrl c stops the code
%; at end prevent console to show results
clc %clear command window
clear %clear workspace window
%press tab to show auto fill code like ctrl space in R
%better to have these at the beginnig of our script
clc;
clear;
close all;
%if you select and press tab it will indent selection
%#ok mean that the error is ok

%right click on command to see help or arguments

%How to run script file line by line?
%first creat a break point (click on lines number on left) then run it then press 'step' icon to run line by line
%if you press 'continue' icon, it will goes to next break point
%if you press 'step in' icon, it goes inside the code of each used function
%in code! cool
%%
%%Session1: Intro

%MATLAB BIBLE: matrix (array) lab. lab is for academic jobs and matrix because all
%data are in matric format
%MATLAB has two face: software and coding
%MATLAB is slower than similar software (C, JAVA, etc.). it is good for academic works
%MATLAB is rapid application development (RAD)
%MATLAB R2013a: R means release in each year they have two releases

%double click on each window maximize and minimize it. in home bar there is
%a layout which can reorganize the windows.
%in current folder window, active directory is shown and matlab will search
%that fist
%workspace window shows variables

A=1
B=[1,2,3]  %[] is a concatnator
%you can omit ,
c=[1,2;3,4] %; means go next line

[A, A] %concating a matrix horizontally
[A; A] %concating a matrix vertically

%if you doing a job takes time it means that you are doing wrong
repmat(D,4,5) %repeating a matrix
A(1,2) %calling i,j=1,2 component
A(2,end) %calling i,j=1,last col component

3:7 %defining a range
3:2:8 %defining a seq by 2 as step length.step in middle

A(2,3:4) %calling a range of component in a matrix
A(1,:) %calling whole of first row equal to 1:end
A(:,1) %calling whole of first col 
A(:,:) %calling all of matrix

A(2) %is second row and first col: is count in col direction

A(:) %changes matrix to a vextor
A(:)' %transposing

sum(1:50) %sum of a range
prod(1:5) %factorial

sum(A) %sum of a matrix cols. toward col direction is default in MATLAB

%defining 3d matrix
c(:,:,1)=[1 2] %first floor
c(:,:,2)=[3 4] %second floor

doc sum %is like ?sum in R for Help
        %write your command but dont close paranthesis and wait, then arg
        %of the command will be shown

sum(A,1) %sum of rows: the result is a row
sum(A,2) %sum of cols: the result is a col
sum(A,3) %sum of floors: the result is a matrix

sum(A(:))%sum of all components

reshape(A,4,3) %change a vector to matrix by 4 row and 3 col
reshape(A,4,[]) %change a vector to matrix by 4 row and any col

diag(A) %shows the diameter elements
diag(diag(A)) %makes all elements of A as zero except diameter

randi([20 100],6,6) %generates random matrix between 20 to 100 as 6*6

%gives interger of 1.3
floor(1.3) 
ceil(1.3) 
round(1.3)

trace(A) %sum of A diameter

det(A) %determinum A

eig(A) %eigen values of A
prod(eig(A))%is same as det(A) 

Filename=num2str(zeros(2,3)) %produce a 2*3 matrix of zeros as str

A=zeros(5,5,2) %generating zero matrix as 5*5*2
A=zeros([5 5 2]) %generating zero matrix as 5*5*2
A=zeros(5) %generating zero matrix as 5*5

ones(5,5,2) %generating 1 matrix as 5*5*2
ones(5)     %generating 1 matrix as 5*5

repmat(0.2,5,5) %generating 0.2 matrix as 5*5
repmat(A,5,5)   %generating A array as a matrix as 5*5

j %reserved for complex number
i %reserved for complex number
i*i %is -1

eye(10) %using same pronounciation to generate I matrix (unit matrix) 
eye(4,7) %generating 1 is dia as possible and others zero

inf %reserved for very big number
inf(3,5) %repmat(inf,3,5)
inf==inf %true
isinf(inf) %true

nan %0/0 which is not a number
nan==nan %is even false :)
isnan(nan) %this is true
%%
%%Session 2: Data types
%check logically: T or F
A==B
A~=B %~ is Contradictory
A<=B
A>=B

%conditions
A>=3 & A<=8 %and
A>=3 | A<=8 %or

%De Morgan's laws: p and q are logical variables
%~(p & q)= ~p | ~q
%~(p | q)= ~p & ~q

find(A<5) %gives the index of elemnts when the condition is true logicals
A(A<5 & A>0) %gives the value of elemnts when the condition is true logicals

mod(A,2) %gives the mod of A devided to 2
sort(A) %sorts the elements of A
class(B) %determine a class of a variable

%data types in MATLAB
%double: Q numbers (real numbers) with high decimal degrees
%single: like doubel with less decimal degrees
realmax %reserved for max double in MATLAB, higher is inf
realmin %reserved for min double in MATLAB, lowere in 0

%charachter:
A='hello man'
double(A) %change char to double
char(b)   %change double to char 
lower(A)  %change all letter to lowercase
upper(A)  %change all letter to uppercase
strcmp(A,B) %check A and B as char as same or not A==B doesnt work on chars

%u+int=unit: unsigned integer:
%if we have n bit we can have 2^n states so:
%unit08: covers 0 to 2^08-1=>0-255 each color a number 0-255
%unit16: covers 0 to 2^16-1=>0-6553
%unit32: covers 0 to 2^32-1
%unit64: covers 0 to 2^64-1
%int08:  covers -2^(08-1) to 2^(08-1)-1=>-128-127
intmax %highest integer
intmin %lowest integer
eps    %reserved for epsilum
imshow(uint8(randi([0 255],1000,1000)));  %showing a random picture

%cell array: having anny kind of data even a matrix 
A={'Ali','Mohammad',23,matrix}
A(1) % is 'Ali'
A{1} %gives content of a cell

%structure data type: a complex type of data having different laye of
%information
%Building structure data type
%x and y are features of s as observations
s.x=[] 
s.y=[]
A=repmat(s,100,1) %A is a strusture of 100s observations which have x and y features
A(1) %gives first row (observation) of A 
%Assigning values to structure data type
A(10).x=10 
A(10).y=5
%%
%Session 3: Control structures in MATLAB
%for loop: appropriate when max of iteration is determined e.g. 50 times
F=zeros(1,10000) %define empty variable before for to provide allocation for saving
%variables as results. it helps run time efficiency.1 here makes vector of
%0s
for i=[1 4 8 3 0]
    disp([i i^2 i^3])
end

%while loop: the max of iteration is not clear
i=1 %first defien i otherwise it condider (-1)^0.5 as predefined i
while(i<5) %inside () is a condition
disp(i^2)
i=i+1; %remember to define i increase
end
%if condition
for i=1:10 %condition should come immidiately after if
    if i==5
        continue   %continue is good for jumping a step
    end
    disp(i)
end

for i=1:10
    if i==5
        break   %terminates loop
    end
    disp(i)
end

%plot: select your variable and go to plot tab and choose desired plot type
A=[3 4 6 7 7]
B=[5 7 8 9 3]
A/B  %/ is matrix division
A./B %./ is for element by element division: A and B should have same length
A*B  %same logic
A.*B %same logic

all(A==B) %if all element are same give 1
any(A==B) %if even one element in each are same gives 1

%run time
tic; %counting sec 
%some coded
toc; %from tic to here 

A(5:10)=[] %elemnt 5th to end will be omitted

disp('This is hello') %showing in the beginning of a code after run
w=input('Enter your age: '); %defining input of a code
disp(['your a is=' num2str(A)]) %good to be used for showing the result of a code
%[] is for pasting str and A as a result of our code
num2str(A) %change num to str and should be used in disp because inside [] the elements should be same
disp({'your a is=' A}) %can be used instead because {} is array and inside it we can have str and double

w=input('Enter your age: ','s'); %means input is a string

mod(n,i) %give mod of n/i

%function is a script which takes input
%name=IsPrime***should be saved in a SEPERATE MFILE with SAME NAME OF FUNCTION
%input=n 
%output=b
function b=IsPrime(n) 
b=true;
    for i=2:floor(sqrt(n))
            if mod(n,i)==0
                b=false;
                break;
            end
    end
end
%we can check our function: 
%it will work ONLY when MATLAB current folder windows contain the function mfile 
IsPrime(4)

%also, we call our function in our code: just becareful about name of file
%and mfile name
m=input("enter input number:")
if IsPrime(m) %the point is that you dont have to consider input only as n it can be m
        disp("input is prime")
        else
        disp("input is not prime")
end

A=[1 3]
i=[5]
B=[A i] %B covers [1 3 5]. this feature is good for loop to add specific results to past defined vector
%%
%Session 4: Functions

%function apply 
f = @(x) x(1,:).^2+x(2,:).^2; %f has two inputs. 
x=[1 2;                   
    2, 3];

output = f(x)  %x rows are inputs so it iterate over x row (first dimesion) itself

%In Cir function we have one input and two output so we can call it:
function [p, a]=Cir(r)
    p=2*pi*r;
    a=pi*r^2;
end
Cir(5) %it shows only first argument here p
[p, a]=Cir(5) %it shows both argument 
[~, a]=Cir(5) %it shows only second argument here a
[~, ~]=Cir(5) %it shows none

%In MatMy we have two inputs and one output
function A=MatMy(m,n)
    if  nargin<1 %nargin means the number of entered args by user
        error('You must enter one arg at least');
    end
    
    if  nargin<2 %nargin means the number of entered args by user
        n=m;
    end
     A=zeros(m,n);
     for i=1:m
         for j=1:n
             A(i,j)=max(i,j);
         end
     end       
end

MatMy(2,3)

%In Prod function we have infinitive inputs to be multiplied to eachother
function A=Prod(varargin) %use varargin for inf inputs: varargin is a cell
    A=1;
    for i=1:numel(varargin)
       A=A.*varargin{i};
    end
end

%In Prodp function we have infinitive inputs would be multiplied to each
%other then power p. you should give p before varargin
function A=Prodp(p,varargin) %use varargin for inf inputs: varargin is a cell
    A=1;
    x=cell2mat(varargin)  %changing cell to matrix
    for i=1:numel(varargin)  %numel give the length of a cell
       A=A*x(i);
    end
    A=A^p;
end

%try to write your code MODULAR: each task done by a specific function
%each function does an specific task
%this approach is very important. first write your script then subdivide it
%into functions

exist('a','var') %we ask q from matlab does any a exist as var? 

%function handles: 3 functionality: 
%1. inherit from other functions and 
fun = @sin;
fminbnd(fun, -3, 3)

%2.make it more easy to call a function
functionX2 = @(x, y) x*2 + y.^2; %it can be run in script

function output = functionX2newface(x, y) %it will run by calling grom another m file
         output = x*2 + y.^2;
end

%3. combination of two above:
ff = @(x) x.^2 + 1;
f_int = integral(ff, -1, 1)
%%
%Session 5: Plotting1
x=1:10;
y=1:10;
z=11:20;
plot(x,y,'o:','Color',[255 0 255]/255,'LineWidth',2,'MarkerSize',3) %desired [R G B] color (balck is [0 0 0] and white is [1 1 1]) by point format: more detail on help of plot
%in plot first you write argument then the value of the argument!
grid on  %showing grid
hold on  %no you can add another plot or even add desired point 
plot(x,z,'r') %with red color
plot(5,5,'s','MarkerSize',10) %add a point to plot
axis equal %fix the axis equally
legend('line1','line2','point')  %having a legend in plot

%how to devide layer of plots 
subplot(2,2,1) %211 means devide layer to 2 rows and one col, so 4 sections. Then, put fig in section 1
plot(x,z,'r') 

subplot(2,2,2) %... put fig in section 2
plot(x,y,'r') 

subplot(2,2,3) %... put fig in section 3
plot(x,z,'r') 

subplot(2,2,4) %... put fig in section 4
plot(x,y,'r') 

col=hsv(5) %produce 5 color. its good for loops!
x=0:0.01:2;
i=0
for a=1:5
    i=i+1; %we need two variable to be used in each loop step so we introduce i
    %i would be used to call hsv for color argument of plot
    y=x.^a; %usig dot in MATLAB: deals with variable arguments one by one istead of dealing as matrix
    plot(x,y,'DisplayName',['a=',num2str(a)],'Color',col(i,:))
    hold on
end
legend show

%how to save plots?

%Go to edit menu>copy figure>paste
%in word: it will have high resolution 
%in Visio, it is editible

%Or save as menu:
%.emf can save fig as a vector readable by word
%.fig only in MATLAB is readable and you can load it without need to the
%source code!!! just double click on it. 

%it plot window, on the most top right (next to cursor icon) there is
%a icon which make it possible to edit plots easily.

%it plot window, on top, there is a insert legend icon
%%
%Session 6: : Plotting2
t=linspace(0,1,1000);  %it generates 1000 num between 0 to 1
%t here is a vectro 1*1000 so t*t is not possible t.*t is possible

%if you select vars in workspace win, and go to plot then you can plot them
%by desired format. In most top left you can chang x-axis and y-axis
%variables as well.

%plot(x,y): in axis all intervals are divided equally
%semilogx(x,y): we zoom in lowest interval in x axis. 0-10 is same as 10-100
%semilogy(x,y): we zoom in lowest interval in y axis. 0-10 is same as 10-100
%loglog(x,y): we zoom in lowest interval in both aixses

y=exp(10*t)
semilogy(t,y) %is the best option when y is exponential
semilogx(t,y) %is worst that plot because you zoon in x axis but y is exponential

x=randn(2,3) %generates 100 rand numbers as a 2*3 matrix
x=randn(1,100) %generates 100 rand numbers as a 1*100 matrix (=vector)

pie(1:3)
colormap gray %changing color of sections of the pie
colormap hsv

explode=[0 0 1 0] %determines which section be departed
pie(1:4,explode)

hist(randn(1,200000),1000) %ploting hist of randn(1,200000) by 1000 bars
histfit(randn(1,200000),1000) %ploting hist and fit norm dist

x=abs(randn(1,200000));
histfit(x,1000,'gamma') %ploting hist and fit gamma dist

hist(x,[1 0 2]) %defining center of a bins

[A C]=hist(x) %takes outputs from hist function instead of plot. A is frequency of each class and C is center of them

%3D plotting: assume you want to plot z=x+3Y
x=[1 2]
y=[4 5 6]
[X Y]=meshgrid(x,y) %first we should consider all possible combinations of x and y in z
                    %so we have 6 combination of x and y. meshgrid generate
                    %them
%here we generate longer x and y:
x=linspace(0,1,20);
y=linspace(0,2,30);
[X, Y]=meshgrid(x,y); %600 combination
Z=X+3*Y;

figure;  %it helpe you to have seperate windows for each plot
surf(X,Y,Z);
colormap hot %it based on z (third dimension value). hotter means higher Z
shading('interp')

subplot(2,2,1);
mesh(X,Y,Z);
colormap hsv
shading('flat')

subplot(2,2,3);
surf(X,Y,Z);
colormap hot 
shading('interp')

subplot(2,2,2);
contour(X,Y,Z);
colormap hsv

subplot(2,2,4);
surfc(X,Y,Z);
colormap hsv

%easy plotting: cool:)
%2D
ezplot('x^2')

%3D
ezsurf('z=x+2*y') %is wrong so:
ezsurf('x+2*y')
%%
%Session 7: Files management1
doc Data and File Management
dir %give current directory features

L=dir
L.name %gives the name of files in directory
{L.name} %it is nicer approach

L=dir('*.m') %it takes only .m files
L=dir('C*')  %it takes only file beginning by C

%its like command prompt of win and linux
cd \  %current folder containg drive
cd .. %current folder dir goes to one level upper
cd 'New Fol' %goes to 'New Fol' excisting folder in current dir

mkdir 'New Fol2' %creates New Fol2 in current dir
rmdir 'New Fol2' %Delet existing New Fol2 in current dir

copyfile('New Fol','Copy New Fol')

isdir('New Fol') %it is folder so it is directory
isdir('Cir.m')   %it is excisting m file so it is not dir   

copyfile('New Fol','Moved Copy New Fol')
movefile('Moved Copy New Fol','New Fol')

%file name construction
matlabroot %the dir of installed matlab
cd(matlabroot) %changed current directory to matlab root dir

toolboxdir('nnet') %toolbox dir

tempdir  %temp file directory
tempname %creates a random name which is unique

filesep %in win filesep is \ but in other sys is / so you can change it

A=[cd filesep 'myfile.m'] %it creates desired directory

[dir name ext]=fileparts(A) %it gives the different elements of a file

uigetfile %easily you can select a file
[A B]=uigetfile %if you select you file, it gives its directory and name 
filename=[B A] %concatinating file name and directory or we can use below command:
fullfile(B,A)

uiputfile %easily you can put a new file

%reading existing excell file xls in current directory: you can easily use
%import icon in above ribbon or:
uiimport %or:
A=xlsread('my excel') 
A=xlsread('my excel','sheet2') %reading desired sheet
A=xlsread('my excel','sheet2','A1:b2') %reading desired range: it is not caps sensitive

[num tx all]=xlsread('my excel','sheet2') %num only reads num in excel, tx only text and all ... :) they are cell type

%writing outpus as excel:
A=[1 2 3]
xlswrite('output',A) %it creats "output" excel file by xls extension
xlswrite('my excel.xlsx',A,'sheet3') %or it saves it in existing "my excel" excel file in new sheet3

%before writing variable from MATLAB workspace we should change in to cell
%1.changing matrix to cell so each element of a matrix is a cell
a=rand(1,10);
b=rand(1,10);
c=rand(1,10);
data=[a; b; c]';
[rowdata coldata]=size(data); %it give raw/col of a data
datanew=mat2cell(data,ones(1,rowdata),ones(1,coldata)); %read help
header={'a' 'b' 'c'};
data_final=[header; datanew]
%2.now exporting is possible
xlswrite('data_final',data_final)
%%
%Session 8: Files management2

%exporting files as .mat files: these files are like excel
a=rand(1,10);
b=rand(1,10);
c=rand(1,10);
data=[a; b; c]';

%you can select variable easily by selecting in workspace and save as .mat
%file by desired name and you can load them:
load('mydata.mat')
load mydata
Q=load('mydata.mat') %Q is a structure

%save files in workspace automatically:
save('mydata_automatic') %by default it save as mat file

Mahdi='is cool'
save('mydata_automatic','Mahdi') %you can define its extension and desired variable 
                                 %note that variable name in in ''
%generating names:
%method1
tempname %creats random name
[dir name]=fileparts(tempname) %it has two outputs 
[~, name]=fileparts(tempname)   %it means only show one output, comma after ~ is importan

%method2
c=clock %it generate time from pc
c(1)    %is year
c(2)    %is month and so on

function name=GetClock
    c=clock;
    Y=c(1);
    M=c(2);
    D=c(3);
   
    h=c(4)
    m=c(5)
    s=floor(c(6))
    
    name=[num2str(Y) '_' num2str(M) '_' num2str(D) '_' ... %... it means that code will continue below
          num2str(h) '_' num2str(m) '_' num2str(s)]
end

%now you can use above function to generate name:
['raster_' GetClock]

%example:
a=1
b=2
c=a+b
filname=['c_' GetClock]
save(filname)

%reading csv(comma seperated value)
csvread
csvwrite
textread
textwrite

A=[1:3; 4:6; 7:9]
dlmwrite('my_datFile.dat',A,'delimiter','\t') %desired delimiter, here is tab for nicer view
DAT_File=dlmread('my_datFile.dat','\t')

read=fileread('my_datFile.dat'); %reading files as str
read(1) %reads firs element of the text file. good for reading text files information

%how to save figure:
My_image=rand(1000);
imshow(rand(1000))
imwrite(My_image,'my_image.png')

%right click on contents in current folder wind and select show in explorer
%then you would be redirected to its containing folder :) cool