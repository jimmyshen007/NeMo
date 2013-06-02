function mail_from_JAKE(status, jobname,msgout)
%mail_from_JAKE sends email on job completion or error interruption

%status='TEST';
%msgout='First attempt.'
%jobname='faux.m'

if status==1
    subject=['[JAKE] Job Completion Notification: ''' jobname ''''];
elseif status==0
    subject=['[JAKE] ERROR: Script ''' jobname ''' aborted'];
else
    subject=status;
end

%Set email
mail='jake.matlab.status@gmail.com';
password='donniedarko2001';

setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props=java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%Send the email

sendmail({'elocastro@gmail.com' '7184064248@tmomail.net'},subject,msgout);
disp('E-mail notification sent.');

end