import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const root=process.argv[2] || path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const out=process.env.ICI_SUBMISSION_ROOT || path.join(root,'submission_final');
const data=JSON.parse(await fs.readFile(path.join(root,'inputs','publication_tables.json'),'utf8'));
const wb=Workbook.create();const map=[];
function letter(n){let s='';while(n>0){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26);}return s;}
for(const spec of data){
 const sh=wb.worksheets.add(spec.name);sh.showGridLines=false;sh.tabColor='#334E68';
 const cols=Math.max(...spec.sections.map(x=>x.headers.length));
 const end=letter(cols);const titleEnd=letter(Math.min(cols,8));
 sh.getRange(`A2:${titleEnd}2`).merge();sh.getRange('A2').values=[[spec.title]];
 sh.getRange(`A2:${titleEnd}2`).format={font:{name:'Arial',size:13,bold:true},rowHeight:26,wrapText:true};
 sh.getRange(`A3:${titleEnd}3`).merge();sh.getRange('A3').values=[[spec.scope]];
 sh.getRange(`A3:${titleEnd}3`).format={font:{name:'Arial',size:10,italic:true},wrapText:true,rowHeight:42};
 let row=5;const widest={};
 for(let si=0;si<spec.sections.length;si++){
  const s=spec.sections[si],c=s.headers.length,last=letter(c),start=row;
  sh.getRange(`A${row}:${last}${row}`).merge();sh.getRange(`A${row}`).values=[[s.title]];
  sh.getRange(`A${row}:${last}${row}`).format={fill:'#E8EDF2',font:{name:'Arial',size:10,bold:true},rowHeight:28,wrapText:true};row++;
  const header=row;sh.getRange(`A${row}:${last}${row}`).values=[s.headers];
  sh.getRange(`A${row}:${last}${row}`).format={fill:'#334E68',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:48,wrapText:true,horizontalAlignment:'center',verticalAlignment:'center'};
  row++;const first=row,finish=row+s.rows.length-1;
  const range=sh.getRange(`A${first}:${last}${finish}`);range.values=s.rows;
  range.format={font:{name:'Arial',size:10,color:'#202020'},rowHeight:24,verticalAlignment:'center',wrapText:true};
  const formats=s.rows.map(r=>r.map((v,j)=>typeof v==='number' ? (Number.isInteger(v)&&(/\(n\)|donors|members|substates|ligands|combinations|sections|cells|spots|Matched library/i.test(s.headers[j]))?'0':(v!==0&&Math.abs(v)<0.001)||/(^P$| q$| P$|GSEA q|Atlas q)/.test(s.headers[j])?'0.000E+00':'0.0000'):'@'));
  range.format.numberFormat=formats;
  for(let j=0;j<c;j++){
   const h=s.headers[j];const textCol=s.rows.some(r=>typeof r[j]==='string'&&!['NA','NE','Yes','No'].includes(r[j]));
   let width=textCol?24:18;
   if(/Gene symbol|Donor$|Donor ID|Ligand|Family$|^Cohort$/.test(h))width=19;
   if(/source|definition|reason|names|members$|subgroup|dimension|scope|Status|Eligibility|Assessment role|Comparison|state|Score representation|Membership comparison/i.test(h)&&textCol)width=38;
   if(h==='Fixed definition')width=64;
   if(h==='Shared members'||h==='Unmeasured fixed members')width=80;
   widest[j+1]=Math.max(widest[j+1]||0,width);
   sh.getRange(`${letter(j+1)}${first}:${letter(j+1)}${finish}`).format.horizontalAlignment=textCol?'left':'center';
   sh.getRange(`${letter(j+1)}${first}:${letter(j+1)}${finish}`).format.borders={left:{style:'hair',color:'#E3E7EB'}};
  }
  // Explicit heights for long genuine table text; no shrinking of body typography.
  for(let i=0;i<s.rows.length;i++){
   let lines=1;
   for(let j=0;j<c;j++)if(typeof s.rows[i][j]==='string'){
    const h=s.headers[j];let chars=30;if(h==='Fixed definition')chars=61;if(h==='Shared members'||h==='Unmeasured fixed members')chars=76;if(h==='Gene symbol'||h==='Donor ID')chars=18;
    lines=Math.max(lines,Math.ceil(s.rows[i][j].length/chars));
   }
   if(lines>1)sh.getRange(`A${first+i}:${last}${first+i}`).format.rowHeight=Math.max(26,lines*13+5);
  }
  row=finish+2;const noteStart=row;
  for(const note of s.notes){sh.getRange(`A${row}:${titleEnd}${row}`).merge();sh.getRange(`A${row}`).values=[[note]];sh.getRange(`A${row}:${titleEnd}${row}`).format={font:{name:'Arial',size:9,color:'#303030'},wrapText:true,rowHeight:Math.max(27,Math.ceil(note.length/160)*13+8),verticalAlignment:'center'};row++;}
  map.push({sheet:spec.name,section:s.title,header,first,last:finish,nrows:s.rows.length,ncols:c,noteStart,noteEnd:row-1});row+=2;
 }
 for(const [j,width] of Object.entries(widest))sh.getRange(`${letter(Number(j))}1:${letter(Number(j))}${row}`).format.columnWidth=width;
 sh.freezePanes.freezeRows(6);
 console.log(spec.name,'rows',row);
}
wb.recalculate();
const errs=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},summary:'Formula error scan'});
await fs.writeFile(path.join(root,'qa','workbook_error_scan.ndjson'),errs.ndjson);
await fs.writeFile(path.join(root,'qa','workbook_cell_ranges.json'),JSON.stringify(map,null,2));
const xlsx=await SpreadsheetFile.exportXlsx(wb);await xlsx.save(path.join(out,'Supplementary_Tables.xlsx'));
// Every worksheet and each distinct section header/footnote receive a visual sample.
for(const spec of data){
 const ranges=map.filter(x=>x.sheet===spec.name);
 for(let j=0;j<ranges.length;j++){
  const m=ranges[j];
  const top=j===0?1:m.header-1;const bottom=Math.min(m.last,m.first+3);
  const preview=await wb.render({sheetName:spec.name,range:`A${top}:${letter(m.ncols)}${bottom}`,scale:1.2,format:'png'});
  await fs.writeFile(path.join(root,'qa',`${spec.name.replace(' ','_')}_${j+1}.png`),new Uint8Array(await preview.arrayBuffer()));
  const notes=await wb.render({sheetName:spec.name,range:`A${m.noteStart}:${letter(Math.min(m.ncols,8))}${m.noteEnd}`,scale:1.2,format:'png'});
  await fs.writeFile(path.join(root,'qa',`${spec.name.replace(' ','_')}_${j+1}_notes.png`),new Uint8Array(await notes.arrayBuffer()));
 }
 console.log('Rendered',spec.name);
}
console.log('Workbook exported and all worksheet sections rendered.');
