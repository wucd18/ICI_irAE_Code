import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
// A legitimate installed module entry must be configured explicitly. No private cache fallback.
if (!process.env.ICI_ARTIFACT_TOOL_MODULE) throw new Error('PORTABILITY_PENDING: ICI_ARTIFACT_TOOL_MODULE is required; see environments/README.md');
const {Workbook,SpreadsheetFile}=await import(pathToFileURL(path.resolve(process.env.ICI_ARTIFACT_TOOL_MODULE)).href);
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const data=JSON.parse(await fs.readFile(path.join(root,'supplement/workbook_input.json'),'utf8'));
const wb=Workbook.create();
const report=[];
for(const d of data){
 const sh=wb.worksheets.add(d.sheet);sh.showGridLines=false;
 const n=d.rows.length,m=d.columns.length;
 sh.getRangeByIndexes(0,0,n+1,m).values=[d.columns,...d.rows];
 const area=sh.getRangeByIndexes(0,0,n+1,m);area.format.font={name:'Arial',size:10};area.format.rowHeight=25;area.format.verticalAlignment='center';
 area.format.borders={preset:'all',style:'thin',color:'#D9D9D9'};
 const head=sh.getRangeByIndexes(0,0,1,m);head.format.fill='#E7EDF2';head.format.font={name:'Arial',size:10,bold:true};head.format.wrapText=true;head.format.rowHeight=44;
 for(let c=0;c<m;c++){
  const vals=d.rows.map(r=>r[c]);const numbers=vals.filter(v=>v!==null).every(v=>typeof v==='number');
  const maxLen=Math.max(d.columns[c].length,...vals.filter(v=>typeof v==='string').map(v=>v.length));
  const col=sh.getRangeByIndexes(0,c,n+1,1);col.format.columnWidthPx=numbers?135:Math.min(340,Math.max(155,maxLen*6));
  if(numbers)sh.getRangeByIndexes(1,c,n,1).setNumberFormat(/(^n_|_n$|size|count)/i.test(d.columns[c])?'0':'0.0000E+00');
  else sh.getRangeByIndexes(1,c,n,1).format.wrapText=true;
 }
 for(let r=0;r<n;r++){
  const height=Math.max(25,...d.rows[r].map((v,c)=>typeof v==='string'?Math.ceil(v.length/(Math.min(340,Math.max(155,Math.max(d.columns[c].length,...d.rows.map(z=>typeof z[c]==='string'?z[c].length:0))*6))/6))*14+8:25));
  sh.getRangeByIndexes(r+1,0,1,m).format.rowHeight=Math.min(130,height);
 }
 sh.freezePanes.freezeRows(1);
 sh.getCell(0,m+1).values=[[d.title]];sh.getCell(0,m+1).format.font={name:'Arial',size:11,bold:true};
 sh.getRangeByIndexes(0,m+1,5,1).format.columnWidthPx=500;
 sh.getCell(2,m+1).values=[[d.source]];sh.getCell(2,m+1).format.wrapText=true;sh.getCell(2,m+1).format.rowHeight=60;
 sh.getCell(4,m+1).values=[[d.note]];sh.getCell(4,m+1).format.wrapText=true;sh.getCell(4,m+1).format.rowHeight=80;
 report.push({sheet:d.sheet,rows:n,columns:m,title:d.title});
}
wb.recalculate();
await fs.mkdir(path.join(root,'review/workbook_previews'),{recursive:true});
for(const d of data){
 const endCol=String.fromCharCode(64+Math.min(5,d.columns.length));
 const img=await wb.render({sheetName:d.sheet,range:`A1:${endCol}${Math.min(d.rows.length+1,7)}`,scale:1.2,format:'png'});
 await fs.writeFile(path.join(root,`review/workbook_previews/${d.sheet}.png`),new Uint8Array(await img.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);await xlsx.save(path.join(root,'supplement/Supplementary_Tables.xlsx'));
const inspect=await wb.inspect({kind:'sheet',include:'id,name',maxChars:2500});
await fs.writeFile(path.join(root,'review/workbook_verification.json'),JSON.stringify({sheets:report,inspect},null,2));
console.log('WORKBOOK_EXPORT_PASS',report.length,'sheets');
