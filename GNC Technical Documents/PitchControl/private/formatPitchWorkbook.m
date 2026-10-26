function formatPitchWorkbook(workbook,sheetNames,sheetTables,formulaCells)
%FORMATPITCHWORKBOOK Native OOXML styling and small executable formula sample.
% Uses MATLAB's bundled Java XML and ZIP support; no Excel/COM dependency.
% Only modifies the newly generated workbook, never source files.
stage = tempname(fileparts(workbook));
mkdir(stage);
cleanup = onCleanup(@() rmdir(stage,'s'));
unzip(workbook,stage);
ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main';
stylesFile = fullfile(stage,'xl','styles.xml');
styles = xmlread(stylesFile);
fonts = styles.getElementsByTagName('fonts').item(0);
fills = styles.getElementsByTagName('fills').item(0);
xfs = styles.getElementsByTagName('cellXfs').item(0);
fontId = fonts.getElementsByTagName('font').getLength;
fillId = fills.getElementsByTagName('fill').getLength;
headerStyle = xfs.getElementsByTagName('xf').getLength;
font = styles.createElementNS(ns,'font');
nameNode = styles.createElementNS(ns,'name'); nameNode.setAttribute('val','Calibri'); font.appendChild(nameNode);
font.appendChild(styles.createElementNS(ns,'b'));
sizeNode = styles.createElementNS(ns,'sz'); sizeNode.setAttribute('val','11'); font.appendChild(sizeNode);
color = styles.createElementNS(ns,'color'); color.setAttribute('rgb','FFFFFFFF'); font.appendChild(color);
fonts.appendChild(font); fonts.setAttribute('count',num2str(fontId+1));
fill = styles.createElementNS(ns,'fill');
pattern = styles.createElementNS(ns,'patternFill'); pattern.setAttribute('patternType','solid');
color = styles.createElementNS(ns,'fgColor'); color.setAttribute('rgb','FF17365D');
pattern.appendChild(color); fill.appendChild(pattern); fills.appendChild(fill);
fills.setAttribute('count',num2str(fillId+1));
xf = styles.createElementNS(ns,'xf');
xf.setAttribute('numFmtId','0'); xf.setAttribute('fontId',num2str(fontId));
xf.setAttribute('fillId',num2str(fillId)); xf.setAttribute('borderId','0');
xf.setAttribute('xfId','0'); xf.setAttribute('applyFont','1'); xf.setAttribute('applyFill','1');
xf.setAttribute('applyAlignment','1');
alignment = styles.createElementNS(ns,'alignment');
alignment.setAttribute('wrapText','1'); alignment.setAttribute('vertical','center');
xf.appendChild(alignment); xfs.appendChild(xf);
bodyStyle = headerStyle+1;
xf = styles.createElementNS(ns,'xf');
xf.setAttribute('numFmtId','0'); xf.setAttribute('fontId','0');
xf.setAttribute('fillId','0'); xf.setAttribute('borderId','0'); xf.setAttribute('xfId','0');
xf.setAttribute('applyAlignment','1');
alignment = styles.createElementNS(ns,'alignment');
alignment.setAttribute('wrapText','1'); alignment.setAttribute('vertical','top');
xf.appendChild(alignment); xfs.appendChild(xf); xfs.setAttribute('count',num2str(bodyStyle+1));
xmlwrite(stylesFile,styles);

for k = 1:numel(sheetNames)
    sheetFile = fullfile(stage,'xl','worksheets',sprintf('sheet%d.xml',k));
    doc = xmlread(sheetFile);
    root = doc.getDocumentElement;
    data = doc.getElementsByTagName('sheetData').item(0);
    oldCols = doc.getElementsByTagName('cols');
    if oldCols.getLength > 0, root.removeChild(oldCols.item(0)); end
    cols = doc.createElementNS(ns,'cols');
    for j = 1:width(sheetTables{k})
        label = sheetTables{k}.Properties.VariableNames{j};
        w = 22;
        if any(strcmp(label,{'TextValue','Source','Details'})), w = 85; end
        if any(strcmp(label,{'Check','Item','Quantity','Confirmation','Status'})), w = 36; end
        col = doc.createElementNS(ns,'col');
        col.setAttribute('min',num2str(j)); col.setAttribute('max',num2str(j));
        col.setAttribute('width',num2str(w)); col.setAttribute('customWidth','1');
        col.setAttribute('style',num2str(bodyStyle)); cols.appendChild(col);
    end
    root.insertBefore(cols,data);
    rows = doc.getElementsByTagName('row');
    firstRow = rows.item(0);
    firstRow.setAttribute('ht','40'); firstRow.setAttribute('customHeight','1');
    cells = firstRow.getElementsByTagName('c');
    for j = 0:cells.getLength-1, cells.item(j).setAttribute('s',num2str(headerStyle)); end
    if any(strcmp(sheetNames{k},{'Inputs_Verification','Quality_Checks','Sensitivity'}))
        for j = 1:rows.getLength-1
            row = rows.item(j); row.setAttribute('ht','46'); row.setAttribute('customHeight','1');
            cells = row.getElementsByTagName('c');
            for m = 0:cells.getLength-1, cells.item(m).setAttribute('s',num2str(bodyStyle)); end
        end
    end
    views = doc.getElementsByTagName('sheetViews');
    if views.getLength == 0
        viewRoot = doc.createElementNS(ns,'sheetViews');
        view = doc.createElementNS(ns,'sheetView'); view.setAttribute('workbookViewId','0');
        viewRoot.appendChild(view);
        formatNodes = doc.getElementsByTagName('sheetFormatPr');
        before = cols;
        if formatNodes.getLength > 0, before = formatNodes.item(0); end
        root.insertBefore(viewRoot,before);
    else
        view = views.item(0).getElementsByTagName('sheetView').item(0);
    end
    view.setAttribute('showGridLines','0');
    pane = doc.createElementNS(ns,'pane');
    pane.setAttribute('ySplit','1'); pane.setAttribute('topLeftCell','A2');
    pane.setAttribute('activePane','bottomLeft'); pane.setAttribute('state','frozen');
    view.insertBefore(pane,view.getFirstChild);
    if doc.getElementsByTagName('autoFilter').getLength == 0
        filter = doc.createElementNS(ns,'autoFilter');
        filter.setAttribute('ref',sprintf('A1:%s%d',excelColumn(width(sheetTables{k})),height(sheetTables{k})+1));
        root.insertBefore(filter,data.getNextSibling);
    end
    if strcmp(sheetNames{k},'Formula_Example')
        cells = data.getElementsByTagName('c');
        for j = 1:size(formulaCells,1)
            found = false;
            for m = 0:cells.getLength-1
                cellNode = cells.item(m);
                if strcmp(char(cellNode.getAttribute('r')),formulaCells{j,1})
                    f = doc.createElementNS(ns,'f');
                    f.appendChild(doc.createTextNode(formulaCells{j,2}));
                    cellNode.insertBefore(f,cellNode.getElementsByTagName('v').item(0));
                    found = true;
                    break
                end
            end
            assert(found,'Pitch:FormulaExport','Formula destination not found: %s.',formulaCells{j,1});
        end
    end
    xmlwrite(sheetFile,doc);
end
bookFile = fullfile(stage,'xl','workbook.xml');
doc = xmlread(bookFile);
calcNodes = doc.getElementsByTagName('calcPr');
if calcNodes.getLength == 0
    calc = doc.createElementNS(ns,'calcPr'); doc.getDocumentElement.appendChild(calc);
else
    calc = calcNodes.item(0);
end
calc.setAttribute('calcMode','auto'); calc.setAttribute('fullCalcOnLoad','1');
xmlwrite(bookFile,doc);
entries = dir(fullfile(stage,'**','*'));
entries = entries(~[entries.isdir]);
relative = cell(numel(entries),1);
for k = 1:numel(entries)
    full = fullfile(entries(k).folder,entries(k).name);
    relative{k} = full(numel(stage)+2:end);
end
archive = fullfile(stage,'formatted.zip');
zip(archive,relative,stage);
movefile(archive,workbook,'f');
end

function label = excelColumn(index)
label = '';
while index > 0
    digit = mod(index-1,26);
    label = [char('A'+digit),label]; %#ok<AGROW>
    index = floor((index-1)/26);
end
end
