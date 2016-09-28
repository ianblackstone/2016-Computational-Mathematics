function [U,G] = bsusurfer(n)
% BSUSURFER Create the adjacency graph of a portion of the websites that
%    have "boisestate" as their domain name.
%    [U,G] = bsusurfer(n) starts at the URL http://www.boisestate.edu and
%    follows Web links until it forms an adjacency graph with n nodes. U =
%    a cell array of n strings, the URLs of the nodes. G = an n-by-n sparse
%    matrix with G(i,j)=1 if node j is linked to node i.
%
%    Example:  [U,G] = bsusurfer(500);
%    See also PAGERANK.
%
%    This function currently has two defects.  (1) The algorithm for
%    finding links is naive.  We just look for the string 'http:'.
%    (2) An attempt to read from a URL that is accessible, but very slow,
%    might take an unacceptably long time to complete.  In some cases,
%    it may be necessary to have the operating system terminate MATLAB.
%    Key words from such URLs can be added to the skip list in bsusurfer.m.
%
%    This is a modified version of the file surfer.m that is provided as
%    part of Cleve Moler's Numerical Computing in MATLAB book.

% Initialize
root = 'http://www.boisestate.edu';
% clf
% shg
% set(gcf,'doublebuffer','on')
% axis([0 n 0 n])
% axis square
% axis ij
% box on
% set(gca,'position',[.12 .20 .78 .78])
% uicontrol('style','frame','units','normal','position',[.01 .09 .98 .07]);
% uicontrol('style','frame','units','normal','position',[.01 .01 .98 .07]);
% t1 = uicontrol('style','text','units','normal','position',[.02 .10 .94 .04], ...
%     'horiz','left');
% t2 = uicontrol('style','text','units','normal','position',[.02 .02 .94 .04], ...
%     'horiz','left');
% slow = uicontrol('style','toggle','units','normal', ...
%     'position',[.01 .24 .07 .05],'string','slow','value',0);
% quit = uicontrol('style','toggle','units','normal', ...
%     'position',[.01 .17 .07 .05],'string','quit','value',0);

U = cell(n,1);
hash = zeros(n,1);
G = logical(sparse(n,n));
m = 1;
U{m} = root;
hash(m) = hashfun(root);

j = 1;
% while j < n & get(quit,'value') == 0
while j < n
    
    % Try to open a page.
    
    if mod(j,10) == 0 || j == 1
        fprintf('number pages processed = %d, %s\n',j,root)
    end
    try
        %set(t1,'string',sprintf('%5d %s',j,U{j}))
        %set(t2,'string','');
        %drawnow
        page = webread(U{j});
    catch
        %set(t1,'string',sprintf('fail: %5d %s',j,U{j}))
        %drawnow
        j = j+1;
        continue
    end
    %if get(slow,'value')
    %   pause(.25)
    %end
    
    % Get the root of the current page
    root = U{j};
    cur_page_root = [];
    re = findstr('/',root);
    if length(re) > 2
        root = root(1:re(3)-1);
        for kk=3:length(re)
            cur_page_root{kk-2} = U{j}(1:re(kk)-1);
        end
        cpn = length(cur_page_root);
    end
    
    % Follow the links from the open page.
    
    %for f = findstr('http:',page);
    href = [findstr('href="',page) findstr('HREF="',page)];
    for f = href
        % A link starts with 'http:' and ends with the next quote.
        
        % e = min([findstr('"',page(f:end)) findstr('''',page(f:end))]);
        % if isempty(e), continue, end
        % url = deblank(page(f:f+e-2));
        % url(url<' ') = '!';   % Nonprintable characters
        % if url(end) == '/', url(end) = []; end
        %
        % % A link starts with 'http:' and ends with the next double quote.
        % f = f+6;
        % e = min(findstr('"',page(f:end)));
        % if isempty(e)
        %    continue
        % end
        % url = deblank(page(f:f+e-2));
        % url(url<' ') = '!';   % Nonprintable characters
        
        % A link starts with 'http:' and ends with the next double quote.
        f = f+6;
        e = min(findstr('"',page(f:end)));
        if isempty(e)
            continue
        end
        url = deblank(page(f:f+e-2));
        url = strtrim(url);      % convert to lowercase
        url(url<' ') = '!';   % Nonprintable characters
        
        if length(url) < 4
            if length(url) == 1 && url(1) == '/'
                url = root;
            else
                continue;
            end
        end
        
        if any(url(1:4) ~= 'http')
            if U{j}(end) ~= '/'
                if ~isempty(findstr(lower(U{j}),'.htm')) || ~isempty(findstr(lower(U{j}),'.html')) || ...
                        ~isempty(findstr(lower(U{j}),'.asp')) || ~isempty(findstr(lower(U{j}),'.shtml')) || ...
                        ~isempty(findstr(lower(U{j}),'.cfm')) || ~isempty(findstr(lower(U{j}),'.shtm')) || ~isempty(findstr(lower(U{j}),'.ico'))
                    if strcmp('../',url(1:3))
                        url = strcat(cur_page_root{max(1,cpn-1)},url(3:end));
                    elseif strcmp('./',url(1:2))
                        url = strcat(cur_page_root{max(1,cpn)},url(2:end));
                    else
                        if url(1) == '/'
                            url = strcat(cur_page_root{cpn},url);
                        else
                            url = strcat(cur_page_root{cpn},'/',url);
                        end
                    end
                else
                    url = strcat(U{j},'/',url);
                end
            else
                if url(1) == '/'
                    url = strcat(U{j},url(2:end));
                else
                    url = strcat(U{j},url);
                end
            end
        end
        
        if url(end) == '/', url(end) = []; end
        
        sw1 = findstr('.boisestate.',lower(url));
        sw2 = findstr('.edu',lower(url));
        if length(sw1) < 1 || length(sw2) < 1
            continue;
        end
        
        % Look for links that should be skipped.
        skips = {'.gif','.jpg','.jpeg','.pdf','.css','.mwc','.ram', ...
            'search.cgi','.cgi','lmscadsi','cybernet','w3.org','google','yahoo', ...
            'scripts','netscape','shockwave','webex','fansonly','www.w3.org',...
            '.ico','..','mailto','.js', '()', 'redirect', 'javascript',...
            '.txt','.wma','.mov','.xml','.doc','.mp3','.xls',...
            '.ppt','.rm','.zip','.asx','.asxp','.ico','.png','.wmv','.nsf','.mpg','/http','/ http',...
            '/ www','.php','shuttle.boisestate',...
            };
        skip = any(url=='!') | any(url=='?') | any(url=='#') | any(url=='$');
        k = 0;
        while ~skip && (k < length(skips))
            k = k+1;
            skip = ~isempty(findstr(lower(url),skips{k}));
        end
        if skip
            %if isempty(findstr(url,'.gif')) & isempty(findstr(url,'.jpg'))
            %   set(t2,'string',sprintf('skip: %s',url))
            %   drawnow
            %   if get(slow,'value')
            %      pause(.25)
            %   end
            %end
            continue
        end
        
        % Remove all double forward slashes, except for the http://
        e = findstr('//',url);
        if numel(e) > 1
            while numel(e) > 1
%                 url = url([1:e(2) min(e(2)+2,length(url)):length(url)]);
                if e(2)+1 == length(url)
                    url = url(1:e(2));
                else
                    url = url([1:e(2) min(e(2)+2,length(url)):length(url)]);
                end
                e = findstr('//',url);
            end
            if url(end) == '/', url(end) = []; end
        end
        
        % Fix the http://boisestate.edu to http://www.boisestate.edu
        e = min(findstr('//',url));
        if numel(e) > 0
            if ( ~isempty(findstr(url(e(1)+2:e(1)+2+10),'boisestate.')) )
                url = [url(1:e(1)+1) 'www.' url(e(1)+2:end)];
            end
        end
        
        % If any webpages are https instead of http, then make them be http
        if isequal(lower(url(1:5)),'https')
            url = url([1:4 6:end]);
        end
%         % It appears that sometimes links start with https when the really
%         % are just http.  Also, some have http when they really should be
%         % https.  I guess browsers deal with this problem separately, but
%         % the matlab code does not.  So we have to manually fix it.
%         if isequal(lower(url(1:5)),'https')
%             % Try to open the page as http.  If this works then alter the
%             % url to just be http.
%             testurl = url([1:4 6:end]);
%             try
%                 testpage = urlread(testurl);
%                 if ~isempty(testpage)
%                     url = testurl;
%                 end
%             catch
%                 j = j+1;
%                 continue;
%             end
%         else
%             % First try to read the page.  If that gives a non-empty
%             % results then everything is assumed to be fine.  Otherwise try
%             % changing it to https.
%             try
%                 testpage = urlread(url);
%                 if isempty(testpage)
%                     url = [url(1:4) 's' url(5:end)];
%                     try
%                         testpage = webread(url);
%                     catch
%                         j = j + 1;
%                         continue;
%                     end
%                 end
%             catch
%                 j = j+1;
%                 continue;
%             end
%         end
                    
        % Check if page is already in url list.
        
        i = 0;
        for k = find(hash(1:m) == hashfun(lower(url)))';
            if isequal(lower(U{k}),lower(url))
                i = k;
                break
            end
        end
        
        % Add a new url to the graph there if are fewer than n.
        
        if (i == 0) && (m < n)
            m = m+1;
            U{m} = url;
            hash(m) = hashfun(lower(url));
            i = m;
        end
        
        % Add a new link.
        
        if i > 0
            G(i,j) = 1;
            %          set(t2,'string',sprintf('%5d %s',i,url))
            %          line(j,i,'marker','.','markersize',6)
            %          drawnow
            %          if get(slow,'value')
            %            pause(.25)
            %          end
            
        end
    end
    
    j = j+1;
end
% delete(t1)
% delete(t2)
% delete(slow)
% set(quit,'string','close','callback','close(gcf)','value',0)



%------------------------

function h = hashfun(url)
% Almost unique numeric hash code for pages already visited.
h = length(url) + 1024*sum(url);
