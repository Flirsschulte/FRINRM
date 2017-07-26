
%**************************************************************************
% Fuzzy Random Impulse Noise Reduction Method 
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
% Inputs:  A = the noisy input image
%          0 = the noise-free image
% Outputs:  Fs = the filtered image 
%**************************************************************************
function Fs = FRINRM(A,O)
   [M,N,DIM] = size(A);
   A = double(A);
   flag = 0;
   if (DIM ~= 1)
       if (sum(sum((A(:,:,1)-A(:,:,2)))) > M) || (sum(sum((A(:,:,1)-A(:,:,3)))) > M)
          Fs = FRINRM_COL(A,O);
          flag = 1;
       else 
          A1 = A(:,:,1);
          clear A;
          A = A1;
       end
   end

   if (flag == 0)
      [M2,N2,DIM2] = size(O);
      if (DIM2 ~= 1)
          O1 = O(:,:,1);
          clear O;
          O = O1;
          clear O1;
      end
      clear O1;

      F = double(A);
      F2 = double(A);
   
      begin = 0;
      loop = 1;
      while (begin == 0)
         F = double(F2);
   
         mem1(M,N,DIM) = 0;
         mem2(M,N,DIM) = 0;

         % neighbourhood size 
         WindowSize  = 1;
         WindowSize2 = 2;
   
         mem1 = Detection1(F,WindowSize);
         mem2 = Detection2(F,WindowSize);
         F2 = Denoise(F,mem1,mem2,WindowSize,WindowSize2);

         if (MSEgray(F,O,3) < MSEgray(F2,O,3)) | (loop == 40)
             begin = 1;
         end
      
         disp(sprintf('==> MSE %2d: %10.5f',loop,MSEgray(F2,O,3)));
         loop = loop +1;
      end
      Fs(:,:,1) = F;
      Fs(:,:,2) = F;
      Fs(:,:,3) = F;
      disp(sprintf('======================='));
      disp(sprintf('Output MSE:  %10.5f',MSEgray(F,O,3)));
      disp(sprintf('Output PSNR: %10.5f',(log(255^2/MSEgray(F,O,3))/log(10)*10)));
      disp(sprintf('======================='));
   end