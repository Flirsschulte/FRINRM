
%**************************************************************************
% Fuzzy Random Impulse Noise Reduction Method for Colour Images 
%   => (Compoent-wise)
%
%  The FRINR method was proposed in: 
%
%  Stefan Schulte, Valérie De Witte, Mike Nachtegael, 
%  Dietrich Van der Weken and  Etienne E. Kerre 
%  Fuzzy Sets and Systems 158(3), 2007 pp. 270-283  
%
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 15/01/06
%
%**************************************************************************
function Fs = FRINRM_COL(A,O)
   [M,N,DIM] = size(A);
   A = double(A);

   A1 = A(:,:,1);
   A2 = A(:,:,2);
   A3 = A(:,:,3);
   
   F1 = double(A1);
   F2 = double(A2);
   F3 = double(A3);
   
   F1b = double(F1);
   F2b = double(F2);
   F3b = double(F3);
   
   begin = 0;
   loop = 1;
   while (begin == 0)
      Fs(:,:,1) = F1b;
      Fs(:,:,2) = F2b;
      Fs(:,:,3) = F3b;
      F1 = double(F1b);
      F2 = double(F2b);
      F3 = double(F3b);
   
      mem1(M,N,DIM) = 0;
      mem2(M,N,DIM) = 0;

      % neighbourhood size 
      WindowSize  = 1;
      WindowSize2 = 2;
   
      mem11 = Detection1(F1,WindowSize);
      mem12 = Detection1(F2,WindowSize);
      mem13 = Detection1(F3,WindowSize);
      
      mem21 = Detection2(F1,WindowSize);
      mem22 = Detection2(F2,WindowSize);
      mem23 = Detection2(F3,WindowSize);
      
      F1b = Denoise(F1,mem11,mem21,WindowSize,WindowSize2);
      F2b = Denoise(F2,mem12,mem22,WindowSize,WindowSize2);
      F3b = Denoise(F3,mem13,mem23,WindowSize,WindowSize2);

      Fs2(:,:,1) = F1b;
      Fs2(:,:,2) = F2b;
      Fs2(:,:,3) = F3b;
      if (MSE(Fs,O,3) < MSE(Fs2,O,3)) | (loop == 40)
          begin = 1;
      end
      
      disp(sprintf('==> MSE %2d: %10.5f',loop,MSE(Fs2,O,3)));
      loop = loop +1;
   end
   disp(sprintf('======================='));
   disp(sprintf('Output MSE:  %10.5f',MSE(Fs,O,3)));
   disp(sprintf('Output PSNR: %10.5f',(log(255^2/MSE(Fs,O,3))/log(10)*10)));
   disp(sprintf('======================='));
