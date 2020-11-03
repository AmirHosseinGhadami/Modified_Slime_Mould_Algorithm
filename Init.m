function Positions = Init(Search_agents_no,dim,ub,lb)

Boundaries_no = size(ub,2);

if Boundaries_no == 1
    Positions = rand(Search_agents_no,dim) .* (ub-lb)+lb;
    
else
    for i=1:dim
        ub_w = ub(i);
        lb_w = lb(i);
        Positions(:,i) = rand(Search_agents_no,1).*(ub_w-lb_w)+lb_w;
    end    
        
end

