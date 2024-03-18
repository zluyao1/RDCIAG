function w = proj2_L1ball0(v,z)
% 函数的功能： 将高维空间点v投影到||w||_1 \leq z的ball中，即在ball中找到一个点与点v之间的欧式距离最小。
% input:
%       v:待投影的矢量(点）
%       z:L1 ball的约束
% output:
%       w:投影后的矢量（点）
% 说明：其中点或者矢量都以列矢量的形式出现
n = length(v);
if norm(v,1)<=z
    w=v;
else
    u=abs(v);
    su=sort(u,1,'descend');
    j=1;
    while  true
        if su(j)-(sum(su(1:j),1)-z)/j<0
            j=j-1;
            break;
        end
        if j==n
            break;
        end
        j=j+1;
    end
    theta=1/j*(sum(su(1:j),1)-z);
    w=bsxfun(@times,sign(v),max(u-theta,0));
end
end
