function [bestParams] = fitLineErrorxy(relX,relY,nguess,minguess,maxguess,tolx,tolfun)
%This function fits a line to x (relX) and y (relY) data by minimizing the
%summed squred distances between the line and the data points (error in 
%both x and y, i.e. not just minimizing deviation in y which is the standard)
    funDiff = @(x) sum(((-x(1)*relX+relY-x(2))/(sqrt((-x(1))^2+1^2))).^2); %Sum of squared distances between a point (normX,normY) and a lines with slope x(1) and intercept x(2)
    
    lb=[-inf -inf];
    ub=[inf inf];
    polySimple=polyfit(relX,relY,1);
    
    guessMin=[minguess*polySimple(1) minguess*polySimple(2)];
    guessMax=[maxguess*polySimple(1) maxguess*polySimple(2)];
    
    options=optimoptions('fmincon','Display','off','TolX',tolx,'TolFun',tolfun);
    numGuesses=nguess;
    paramsAll=zeros(numGuesses,numel(lb));
    resnormAll=zeros(numGuesses,1);
    for j= 1:numGuesses
        
        params0=(guessMin+rand(1,numel(lb)).*(guessMax-guessMin));    
        [paramsCurr,resnorm] = fmincon(funDiff,params0,[],[],[],[],lb,ub,[],options);
        resnormAll(j)=resnorm;
        paramsAll(j,:)=paramsCurr;
        
    end
    
    [~,bestind]=min(resnormAll);
    bestParams=paramsAll(bestind,:);

end

