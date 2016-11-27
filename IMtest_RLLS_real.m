% Information test for RL
clear all global;
clear all;
importfile('./WSRL.mat');
%  importfile('WSRL.mat');
globalVar;

notifyMail('set','amyeuphich@gmail.com','sntal2908');
isFixedUturn = false;
value = IMtest(Op.x);
result =  sprintf(' Information matrix equality test (RL-LS-real): %f\n', value);
try
   notifyMail('send', result);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end