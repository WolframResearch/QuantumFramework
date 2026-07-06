import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(__dirname, "..");
const DATA_DIR = path.join(ROOT, "data");
const OUTPUT_DIR = path.join(ROOT, "outputs", "quantum_textbook_survey_20260705");
const OUTPUT_FILE = path.join(OUTPUT_DIR, "quantum_textbook_survey_public_evidence.xlsx");

async function readData() {
  const read = async (name) => JSON.parse(await fs.readFile(path.join(DATA_DIR, name), "utf8"));
  return {
    universities: await read("universities_arwu_2025_top100.json"),
    courses: await read("courses.json"),
    textbooks: await read("textbooks.json"),
    canonicalBooks: await read("canonical_books.json"),
    gaps: await read("gaps.json"),
    candidates: await read("search_candidates.json"),
    summary: await read("summary.json"),
  };
}

function colLetter(index) {
  let n = index + 1;
  let s = "";
  while (n > 0) {
    const m = (n - 1) % 26;
    s = String.fromCharCode(65 + m) + s;
    n = Math.floor((n - m) / 26);
  }
  return s;
}

function sheetRange(rowCount, colCount) {
  return `A1:${colLetter(colCount - 1)}${Math.max(1, rowCount)}`;
}

function rowsFromObjects(rows, headers) {
  return rows.map((row) => headers.map((header) => row[header] ?? ""));
}

function setHeaderStyle(range) {
  range.format.fill.color = "#1F4E79";
  range.format.font.color = "#FFFFFF";
  range.format.font.bold = true;
  range.format.horizontalAlignment = "center";
  range.format.verticalAlignment = "center";
  range.format.wrapText = true;
}

function setBodyStyle(range) {
  range.format.font.name = "Aptos";
  range.format.font.size = 10;
  range.format.verticalAlignment = "top";
  range.format.wrapText = true;
}

function addDataSheet(workbook, name, headers, rows, widths = {}) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const values = [headers, ...rowsFromObjects(rows, headers)];
  const range = sheet.getRange(sheetRange(values.length, headers.length));
  range.values = values;
  setBodyStyle(range);
  setHeaderStyle(sheet.getRange(`A1:${colLetter(headers.length - 1)}1`));
  sheet.getRange(sheetRange(values.length, headers.length)).format.borders = {
    insideHorizontal: { style: "thin", color: "#E6EAF0" },
    top: { style: "thin", color: "#AAB7C4" },
    bottom: { style: "thin", color: "#AAB7C4" },
  };
  sheet.freezePanes.freezeRows(1);
  for (let c = 0; c < headers.length; c += 1) {
    const width = widths[headers[c]] ?? 18;
    sheet.getRangeByIndexes(0, c, values.length, 1).format.columnWidth = width;
  }
  sheet.getRange(`A1:${colLetter(headers.length - 1)}1`).format.rowHeight = 34;
  if (values.length > 1) {
    sheet.getRange(`A2:${colLetter(headers.length - 1)}${values.length}`).format.rowHeight = 42;
  }
  return sheet;
}

function addReadme(workbook, data) {
  const sheet = workbook.worksheets.add("README");
  sheet.showGridLines = false;
  const generatedAt = data.summary.generated_at.replace("T", " ").replace("Z", " UTC");
  const rows = [
    ["Quantum Textbook Survey - Public Evidence Snapshot", ""],
    ["Generated", generatedAt],
    ["Ranking seed", "ARWU 2025 top 100 from ShanghaiRanking public payload"],
    ["Ranking URL", data.summary.ranking_url],
    ["Scope", "BSc/MSc-level quantum-related course textbook evidence found in public web sources."],
    ["Important limitation", "This is not a complete global census. It is a traceable public-evidence register; LMS, bookstore, and departmental private data remain unavailable without access or outreach."],
    ["Search limitation", "Batch search access was rate-limited. The workbook therefore separates verified rows from candidates and records gaps for unverified universities."],
    ["Confidence A", "Official current/course syllabus or course page explicitly matched a textbook or reading."],
    ["Confidence B", "Official academic source matched a course/readings page but with weaker textbook wording."],
    ["Confidence C", "Academic domain or official-looking page with textbook evidence from snippet/source metadata."],
    ["Confidence D", "Third-party or archived/candidate source with textbook evidence."],
    ["Confidence E", "Low-confidence or non-course source; do not use in headline counts."],
    ["Recommended next step", "For universities in Gaps, check official reading-list systems, bookstores by course code, department course pages, and then email departments/instructors."],
  ];
  const range = sheet.getRange(`A1:B${rows.length}`);
  range.values = rows;
  setBodyStyle(range);
  sheet.getRange("A1:B1").merge();
  sheet.getRange("A1").format.fill.color = "#16324F";
  sheet.getRange("A1").format.font.color = "#FFFFFF";
  sheet.getRange("A1").format.font.bold = true;
  sheet.getRange("A1").format.font.size = 16;
  sheet.getRange("A1").format.rowHeight = 32;
  sheet.getRange(`A2:A${rows.length}`).format.font.bold = true;
  sheet.getRange(`A2:A${rows.length}`).format.fill.color = "#EAF0F6";
  sheet.getRange(`A1:B${rows.length}`).format.borders = { preset: "all", style: "thin", color: "#D6DEE8" };
  sheet.getRange("A:A").format.columnWidth = 24;
  sheet.getRange("B:B").format.columnWidth = 120;
  sheet.freezePanes.freezeRows(1);
}

function addSummary(workbook, data) {
  const sheet = workbook.worksheets.add("Summary");
  sheet.showGridLines = false;
  const metrics = [
    ["Metric", "Value"],
    ["Universities in ranking seed", data.summary.university_count],
    ["Search candidates retained", data.summary.candidate_count],
    ["Course candidates checked", data.summary.course_candidate_count],
    ["Textbook matches, all confidence levels", data.summary.textbook_match_count],
    ["Verified textbook matches, confidence A-C", data.summary.verified_textbook_match_count],
    ["Universities with verified matches", data.summary.universities_with_verified_matches],
    ["Universities requiring manual follow-up", data.summary.universities_with_gaps],
  ];
  sheet.getRange(`A1:B${metrics.length}`).values = metrics;
  setBodyStyle(sheet.getRange(`A1:B${metrics.length}`));
  setHeaderStyle(sheet.getRange("A1:B1"));
  sheet.getRange(`A1:B${metrics.length}`).format.borders = { preset: "all", style: "thin", color: "#D6DEE8" };
  sheet.getRange("A:A").format.columnWidth = 44;
  sheet.getRange("B:B").format.columnWidth = 18;

  const topHeaders = [
    "canonical_title",
    "canonical_authors",
    "category",
    "used_by_university_count",
    "used_in_course_count",
    "confidence_a_to_c_count",
  ];
  const topRows = [topHeaders, ...rowsFromObjects(data.canonicalBooks.slice(0, 15), topHeaders)];
  sheet.getRange(`D1:${colLetter(3 + topHeaders.length - 1)}${topRows.length}`).values = topRows;
  setBodyStyle(sheet.getRange(`D1:${colLetter(3 + topHeaders.length - 1)}${topRows.length}`));
  setHeaderStyle(sheet.getRange(`D1:${colLetter(3 + topHeaders.length - 1)}1`));
  sheet.getRange(`D1:${colLetter(3 + topHeaders.length - 1)}${topRows.length}`).format.borders = { preset: "all", style: "thin", color: "#D6DEE8" };
  const summaryWidths = [34, 34, 22, 18, 18, 18];
  for (let i = 0; i < topHeaders.length; i += 1) {
    sheet.getRangeByIndexes(0, 3 + i, topRows.length, 1).format.columnWidth = summaryWidths[i];
  }

  const noteRows = [
    ["Interpretation"],
    ["Use confidence A-C for reliable counts. D/E rows are retained only as leads for manual verification."],
    ["Because most universities hide syllabi or reading lists behind LMS/bookstore systems, the Gaps sheet is part of the result, not a failure state."],
  ];
  sheet.getRange("A11:H13").values = [
    [noteRows[0][0], "", "", "", "", "", "", ""],
    [noteRows[1][0], "", "", "", "", "", "", ""],
    [noteRows[2][0], "", "", "", "", "", "", ""],
  ];
  sheet.getRange("A11:H11").merge();
  sheet.getRange("A12:H12").merge();
  sheet.getRange("A13:H13").merge();
  sheet.getRange("A11").format.fill.color = "#F3B23C";
  sheet.getRange("A11").format.font.bold = true;
  sheet.getRange("A12:A13").format.fill.color = "#FFF7E6";
  sheet.getRange("A12:A13").format.wrapText = true;
  sheet.freezePanes.freezeRows(1);
}

async function main() {
  const data = await readData();
  const workbook = Workbook.create();

  addReadme(workbook, data);
  addSummary(workbook, data);

  addDataSheet(
    workbook,
    "Universities",
    ["rank", "university", "country", "regional_rank", "ranking_source", "ranking_year", "ranking_url", "arwu_profile_url", "overall_score"],
    data.universities,
    {
      rank: 10,
      university: 42,
      country: 24,
      ranking_url: 58,
      arwu_profile_url: 58,
    },
  );

  addDataSheet(
    workbook,
    "Courses",
    ["university", "rank", "country", "course_code", "course_title", "level", "degree_context", "department", "academic_year", "term", "source_url", "source_type", "source_status", "confidence", "notes"],
    data.courses,
    {
      university: 34,
      course_title: 48,
      level: 24,
      source_url: 72,
      source_status: 34,
      notes: 30,
    },
  );

  addDataSheet(
    workbook,
    "Textbooks",
    ["university", "rank", "country", "course_code", "course_title", "level", "textbook_title", "authors", "edition", "publisher", "year", "isbn", "required_or_recommended", "source_evidence_note", "source_url", "confidence", "category", "matched_terms"],
    data.textbooks,
    {
      university: 34,
      course_title: 46,
      textbook_title: 44,
      authors: 38,
      required_or_recommended: 28,
      source_evidence_note: 48,
      source_url: 72,
      matched_terms: 34,
    },
  );

  addDataSheet(
    workbook,
    "Canonical_Books",
    ["canonical_title", "canonical_authors", "category", "used_by_university_count", "used_in_course_count", "bsc_count", "msc_count", "unknown_level_count", "confidence_a_to_c_count", "source_url_count"],
    data.canonicalBooks,
    {
      canonical_title: 44,
      canonical_authors: 42,
      category: 26,
    },
  );

  addDataSheet(
    workbook,
    "Gaps",
    ["university", "rank", "country", "course_or_program", "missing_item", "reason", "attempted_sources", "recommended_next_action"],
    data.gaps,
    {
      university: 36,
      course_or_program: 28,
      missing_item: 34,
      reason: 56,
      recommended_next_action: 90,
    },
  );

  addDataSheet(
    workbook,
    "Search_Candidates",
    ["university", "rank", "country", "search_rank", "score", "title", "url", "host", "snippet", "query"],
    data.candidates,
    {
      university: 34,
      title: 60,
      url: 76,
      snippet: 90,
      query: 54,
    },
  );

  const errorScan = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 300 },
    summary: "final formula error scan",
  });
  console.log(errorScan.ndjson);

  await fs.mkdir(OUTPUT_DIR, { recursive: true });
  for (const sheetName of ["README", "Summary", "Universities", "Courses", "Textbooks", "Canonical_Books", "Gaps", "Search_Candidates"]) {
    const preview = await workbook.render({ sheetName, range: "A1:H30", scale: 1, format: "png" });
    const bytes = new Uint8Array(await preview.arrayBuffer());
    await fs.writeFile(path.join(OUTPUT_DIR, `${sheetName}.png`), bytes);
  }

  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(OUTPUT_FILE);
  console.log(OUTPUT_FILE);
}

main().catch((error) => {
  console.error(error.stack || error.message);
  process.exit(1);
});
