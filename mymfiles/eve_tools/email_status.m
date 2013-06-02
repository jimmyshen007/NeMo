function email_status(status, jobname, msgout, recip)
%email_status sends email on job completion or error interruption, or
%whatever other occasion is warranted by 'status'; recip is input as string
%cell array

eve_config;

% Auto-recipients set here if nargin < 4

if nargin<4
    recip(1)='elocastro@gmail.com';
    recip(2)='7184064248@tmomail.net';
end
    
S=size(recip);
if S(1) > S(2)
    R=S(1);
else
    R=S(2);
end

if status==1
    subject=['[' hostname '] Job Completion Notification: ''' jobname ''''];
elseif status==0
    subject=['[' hostname '] ERROR: Script ''' jobname ''' aborted'];
else
    subject=int2str(status);
end


%Send the email

for r=1:R
    sendmail(recip{r},subject,msgout);
end
disp('E-mail notification sent.');

end