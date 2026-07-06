import fs from "node:fs/promises";
import path from "node:path";
import vm from "node:vm";
import crypto from "node:crypto";
import http from "node:http";
import https from "node:https";
import { URL } from "node:url";

const ROOT = path.resolve(new URL("../", import.meta.url).pathname);
const DATA_DIR = path.join(ROOT, "data");
const CACHE_DIR = path.join(DATA_DIR, "cache");
const SEARCH_CACHE_DIR = path.join(CACHE_DIR, "search");
const BRAVE_CACHE_DIR = path.join(CACHE_DIR, "brave_search");
const PAGE_CACHE_DIR = path.join(CACHE_DIR, "pages");

const RANKING_PAGE = "https://www.shanghairanking.com/rankings/arwu/2025";
const USER_AGENT =
  "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 " +
  "(KHTML, like Gecko) Chrome/126.0 Safari/537.36";

const NEGATIVE_HOSTS = [
  "wikipedia.org",
  "topuniversities.com",
  "timeshighereducation.com",
  "shanghairanking.com",
  "coursehero.com",
  "chegg.com",
  "quizlet.com",
  "studocu.com",
  "books.google.com",
  "reddit.com",
  "youtube.com",
  "facebook.com",
  "linkedin.com",
  "ratemyprofessors.com",
  "researchgate.net",
];

const BOOKS = [
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "C. Cohen-Tannoudji; B. Diu; F. Laloe",
    category: "Quantum mechanics",
    patterns: [
      /Cohen[-\s]?Tannoudji/i,
      /\bCTDL\b/i,
      /Diu,\s*(?:and\s*)?Lalo/i,
    ],
  },
  {
    canonical_title: "Introduction to Quantum Mechanics",
    canonical_authors: "David J. Griffiths",
    category: "Quantum mechanics",
    patterns: [/Griffiths/i, /Introduction to Quantum Mechanics/i],
  },
  {
    canonical_title: "Principles of Quantum Mechanics",
    canonical_authors: "R. Shankar",
    category: "Quantum mechanics",
    patterns: [/Shankar/i, /Principles of Quantum Mechanics/i],
  },
  {
    canonical_title: "Modern Quantum Mechanics",
    canonical_authors: "J. J. Sakurai; Jim Napolitano",
    category: "Quantum mechanics",
    patterns: [/Sakurai/i, /Modern Quantum Mechanics/i],
  },
  {
    canonical_title: "Quantum Mechanics: Concepts and Applications",
    canonical_authors: "Nouredine Zettili",
    category: "Quantum mechanics",
    patterns: [/Zettili/i, /Concepts and Applications/i],
  },
  {
    canonical_title: "A Modern Approach to Quantum Mechanics",
    canonical_authors: "John S. Townsend",
    category: "Quantum mechanics",
    patterns: [/Townsend/i, /Modern Approach to Quantum Mechanics/i],
  },
  {
    canonical_title: "Quantum Mechanics: A Paradigms Approach",
    canonical_authors: "David H. McIntyre",
    category: "Quantum mechanics",
    patterns: [/McIntyre/i, /Paradigms Approach/i],
  },
  {
    canonical_title: "Quantum Physics",
    canonical_authors: "Stephen Gasiorowicz",
    category: "Quantum mechanics",
    patterns: [/Gasiorowicz/i],
  },
  {
    canonical_title: "Introductory Quantum Mechanics",
    canonical_authors: "Richard L. Liboff",
    category: "Quantum mechanics",
    patterns: [/Liboff/i],
  },
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "B. H. Bransden; C. J. Joachain",
    category: "Quantum mechanics",
    patterns: [/Bransden/i, /Joachain/i],
  },
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "Franz Schwabl",
    category: "Quantum mechanics",
    patterns: [/Schwabl/i],
  },
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "Albert Messiah",
    category: "Quantum mechanics",
    patterns: [/Messiah/i],
  },
  {
    canonical_title: "Quantum Mechanics: A Modern Development",
    canonical_authors: "Leslie E. Ballentine",
    category: "Quantum mechanics",
    patterns: [/Ballentine/i],
  },
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "Eugen Merzbacher",
    category: "Quantum mechanics",
    patterns: [/Merzbacher/i],
  },
  {
    canonical_title: "Quantum Mechanics",
    canonical_authors: "Leonard I. Schiff",
    category: "Quantum mechanics",
    patterns: [/Schiff/i],
  },
  {
    canonical_title: "Quantum Computation and Quantum Information",
    canonical_authors: "Michael A. Nielsen; Isaac L. Chuang",
    category: "Quantum information",
    patterns: [/Nielsen/i, /Chuang/i, /Quantum Computation and Quantum Information/i],
  },
  {
    canonical_title: "Quantum Computer Science",
    canonical_authors: "N. David Mermin",
    category: "Quantum information",
    patterns: [/Mermin/i, /Quantum Computer Science/i],
  },
  {
    canonical_title: "Quantum Computation lecture notes",
    canonical_authors: "John Preskill",
    category: "Quantum information",
    patterns: [/Preskill/i],
  },
  {
    canonical_title: "Quantum Computing: A Gentle Introduction",
    canonical_authors: "Eleanor Rieffel; Wolfgang Polak",
    category: "Quantum information",
    patterns: [/Rieffel/i, /Polak/i, /Gentle Introduction/i],
  },
  {
    canonical_title: "Quantum Information Theory",
    canonical_authors: "Mark M. Wilde",
    category: "Quantum information",
    patterns: [/Wilde/i, /Quantum Information Theory/i],
  },
  {
    canonical_title: "The Theory of Quantum Information",
    canonical_authors: "John Watrous",
    category: "Quantum information",
    patterns: [/Watrous/i, /Theory of Quantum Information/i],
  },
  {
    canonical_title: "Quantum Computing for Computer Scientists",
    canonical_authors: "Noson S. Yanofsky; Mirco A. Mannucci",
    category: "Quantum information",
    patterns: [/Yanofsky/i, /Mannucci/i],
  },
  {
    canonical_title: "Quantum Optics",
    canonical_authors: "Marlan O. Scully; M. Suhail Zubairy",
    category: "Quantum optics",
    patterns: [/Scully/i, /Zubairy/i],
  },
  {
    canonical_title: "Introductory Quantum Optics",
    canonical_authors: "Christopher C. Gerry; Peter L. Knight",
    category: "Quantum optics",
    patterns: [/Gerry/i, /Knight/i, /Introductory Quantum Optics/i],
  },
  {
    canonical_title: "Quantum Optics",
    canonical_authors: "Daniel F. Walls; Gerard J. Milburn",
    category: "Quantum optics",
    patterns: [/Walls/i, /Milburn/i],
  },
  {
    canonical_title: "Optical Coherence and Quantum Optics",
    canonical_authors: "Leonard Mandel; Emil Wolf",
    category: "Quantum optics",
    patterns: [/Mandel/i, /Wolf/i, /Optical Coherence/i],
  },
  {
    canonical_title: "The Quantum Theory of Light",
    canonical_authors: "Rodney Loudon",
    category: "Quantum optics",
    patterns: [/Loudon/i, /Quantum Theory of Light/i],
  },
  {
    canonical_title: "Quantum Theory of Many-Particle Systems",
    canonical_authors: "Alexander L. Fetter; John D. Walecka",
    category: "Many-body quantum theory",
    patterns: [/Fetter/i, /Walecka/i],
  },
  {
    canonical_title: "Condensed Matter Field Theory",
    canonical_authors: "Alexander Altland; Ben Simons",
    category: "Many-body quantum theory",
    patterns: [/Altland/i, /Simons/i, /Condensed Matter Field Theory/i],
  },
  {
    canonical_title: "An Introduction to Quantum Field Theory",
    canonical_authors: "Michael E. Peskin; Daniel V. Schroeder",
    category: "Quantum field theory",
    patterns: [/Peskin/i, /Schroeder/i],
  },
  {
    canonical_title: "Quantum Field Theory in a Nutshell",
    canonical_authors: "A. Zee",
    category: "Quantum field theory",
    patterns: [/\bZee\b/i, /Quantum Field Theory in a Nutshell/i],
  },
  {
    canonical_title: "Quantum Field Theory and the Standard Model",
    canonical_authors: "Matthew D. Schwartz",
    category: "Quantum field theory",
    patterns: [/Schwartz/i, /Quantum Field Theory and the Standard Model/i],
  },
  {
    canonical_title: "The Quantum Theory of Fields",
    canonical_authors: "Steven Weinberg",
    category: "Quantum field theory",
    patterns: [/Weinberg/i, /Quantum Theory of Fields/i],
  },
  {
    canonical_title: "Advanced Quantum Mechanics",
    canonical_authors: "J. J. Sakurai",
    category: "Advanced quantum mechanics",
    patterns: [/Advanced Quantum Mechanics/i],
  },
];

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function sha1(input) {
  return crypto.createHash("sha1").update(input).digest("hex");
}

function decodeEntities(input = "") {
  return input
    .replace(/&#(\d+);/g, (_, n) => String.fromCharCode(Number(n)))
    .replace(/&#x([0-9a-f]+);/gi, (_, n) => String.fromCharCode(parseInt(n, 16)))
    .replace(/&amp;/g, "&")
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&nbsp;/g, " ");
}

function stripTags(input = "") {
  return decodeEntities(input.replace(/<[^>]*>/g, " ")).replace(/\s+/g, " ").trim();
}

function cleanText(input = "") {
  return decodeEntities(input)
    .replace(/<script[\s\S]*?<\/script>/gi, " ")
    .replace(/<style[\s\S]*?<\/style>/gi, " ")
    .replace(/<[^>]+>/g, " ")
    .replace(/\s+/g, " ")
    .trim();
}

function safeCsv(value) {
  if (value == null) return "";
  const s = String(value);
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

async function writeCsv(filePath, rows, headers) {
  const lines = [headers.join(",")];
  for (const row of rows) {
    lines.push(headers.map((header) => safeCsv(row[header])).join(","));
  }
  await fs.writeFile(filePath, lines.join("\n") + "\n", "utf8");
}

async function fetchUrl(url, options = {}) {
  const {
    timeoutMs = 25000,
    maxRedirects = 5,
    maxBytes = 2_500_000,
    cacheDir = null,
    cacheKey = url,
  } = options;

  if (cacheDir) {
    await fs.mkdir(cacheDir, { recursive: true });
    const cachePath = path.join(cacheDir, `${sha1(cacheKey)}.txt`);
    try {
      const cached = await fs.readFile(cachePath, "utf8");
      const parsed = JSON.parse(cached);
      if (parsed.status !== 0) return parsed;
    } catch {}
  }

  async function inner(target, redirectsLeft) {
    return new Promise((resolve, reject) => {
      const parsed = new URL(target);
      const lib = parsed.protocol === "http:" ? http : https;
      const req = lib.request(
        parsed,
        {
          headers: {
            "User-Agent": USER_AGENT,
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Encoding": "identity",
          },
          timeout: timeoutMs,
        },
        (res) => {
          const status = res.statusCode ?? 0;
          const location = res.headers.location;
          if ([301, 302, 303, 307, 308].includes(status) && location && redirectsLeft > 0) {
            res.resume();
            const next = new URL(location, target).href;
            resolve(inner(next, redirectsLeft - 1));
            return;
          }

          const chunks = [];
          let total = 0;
          res.on("data", (chunk) => {
            total += chunk.length;
            if (total <= maxBytes) {
              chunks.push(chunk);
            }
          });
          res.on("end", () => {
            const body = Buffer.concat(chunks).toString("utf8");
            resolve({
              url: target,
              status,
              contentType: res.headers["content-type"] ?? "",
              truncated: total > maxBytes,
              body,
            });
          });
        },
      );
      req.on("timeout", () => {
        req.destroy(new Error(`Timeout fetching ${target}`));
      });
      req.on("error", reject);
      req.end();
    });
  }

  try {
    const result = await inner(url, maxRedirects);
    if (cacheDir) {
      const cachePath = path.join(cacheDir, `${sha1(cacheKey)}.txt`);
      await fs.writeFile(cachePath, JSON.stringify(result), "utf8");
    }
    return result;
  } catch (error) {
    const result = {
      url,
      status: 0,
      contentType: "",
      truncated: false,
      body: "",
      error: error.message,
    };
    if (cacheDir) {
      const cachePath = path.join(cacheDir, `${sha1(cacheKey)}.txt`);
      await fs.writeFile(cachePath, JSON.stringify(result), "utf8");
    }
    return result;
  }
}

async function fetchArwuPayload() {
  const page = await fetchUrl(RANKING_PAGE, {
    cacheDir: CACHE_DIR,
    cacheKey: "arwu-2025-page",
    maxBytes: 1_000_000,
  });
  const payloadMatch = page.body.match(/["']([^"']*_nuxt\/static\/[^"']*\/rankings\/arwu\/2025\/payload\.js)["']/);
  if (!payloadMatch) {
    throw new Error("Could not find ARWU payload.js path in ranking page.");
  }
  const payloadUrl = new URL(payloadMatch[1], RANKING_PAGE).href;
  const payload = await fetchUrl(payloadUrl, {
    cacheDir: CACHE_DIR,
    cacheKey: "arwu-2025-payload",
    maxBytes: 2_000_000,
  });
  let captured;
  vm.runInNewContext(payload.body, {
    __NUXT_JSONP__: (_route, data) => {
      captured = data;
    },
  });
  const data = captured?.data?.[0];
  if (!data?.univList?.length) {
    throw new Error("ARWU payload did not contain univList.");
  }
  return {
    payloadUrl,
    intro: data.introWeb,
    universities: data.univList.slice(0, 100).map((row) => ({
      rank: row.ranking,
      university: row.univNameEn,
      country: row.region,
      regional_rank: row.regionRanking,
      ranking_source: "ARWU",
      ranking_year: 2025,
      ranking_url: RANKING_PAGE,
      arwu_profile_url: `https://www.shanghairanking.com/universities/${row.univUp}`,
      overall_score: row.score,
    })),
  };
}

function decodeDdgUrl(href) {
  const decodedHref = decodeEntities(href);
  const absolute = decodedHref.startsWith("//")
    ? `https:${decodedHref}`
    : new URL(decodedHref, "https://duckduckgo.com").href;
  try {
    const url = new URL(absolute);
    return url.searchParams.get("uddg") || absolute;
  } catch {
    return absolute;
  }
}

function parseDdg(html) {
  const blocks = html.split(/<div class="result /).slice(1);
  const results = [];
  for (const block of blocks) {
    const titleMatch = block.match(/<a[^>]+class="result__a"[^>]+href="([^"]+)"[^>]*>([\s\S]*?)<\/a>/);
    if (!titleMatch) continue;
    const snippetMatch = block.match(/<a[^>]+class="result__snippet"[^>]*>([\s\S]*?)<\/a>/);
    const url = decodeDdgUrl(titleMatch[1]);
    const title = stripTags(titleMatch[2]);
    const snippet = stripTags(snippetMatch?.[1] ?? "");
    results.push({ title, url, snippet });
  }
  return results;
}

function decodeJsString(input = "") {
  try {
    return JSON.parse(`"${input.replace(/\n/g, "\\n")}"`);
  } catch {
    return input
      .replace(/\\"/g, '"')
      .replace(/\\u003C/g, "<")
      .replace(/\\u003E/g, ">")
      .replace(/\\u002F/g, "/");
  }
}

function parseBrave(html) {
  const results = [];
  const seen = new Set();
  const re =
    /(?:^|[,{])title:"((?:\\.|[^"\\])*)",url:"((?:\\.|[^"\\])*)",full_title:(?:void 0|"((?:\\.|[^"\\])*)"),description:"((?:\\.|[^"\\])*)"/g;
  let match;
  while ((match = re.exec(html)) && results.length < 25) {
    const title = stripTags(decodeJsString(match[3] || match[1]));
    const url = decodeJsString(match[2]);
    const snippet = stripTags(decodeJsString(match[4]));
    if (!/^https?:\/\//.test(url) || seen.has(url)) continue;
    seen.add(url);
    results.push({ title, url, snippet });
  }
  return results;
}

function hostOf(url) {
  try {
    return new URL(url).hostname.replace(/^www\./, "");
  } catch {
    return "";
  }
}

function scoreSearchResult(result) {
  const host = hostOf(result.url);
  const text = `${result.title} ${result.url} ${result.snippet}`.toLowerCase();
  let score = 0;
  if (NEGATIVE_HOSTS.some((negative) => host.endsWith(negative))) score -= 20;
  if (/quantum/.test(text)) score += 4;
  if (/syllabus|reading list|readings|textbook|course outline|module|course catalog|opencourseware|lecture notes/.test(text)) score += 6;
  if (/textbook|required text|recommended text|isbn|book/.test(text)) score += 5;
  if (/\.edu$|\.edu\.|\.ac\.uk$|\.edu\.au$|\.edu\.cn$|\.ac\.jp$|\.ac\.kr$|\.edu\.sg$/.test(host)) score += 4;
  if (/pdf$/i.test(result.url)) score += 1;
  if (/coursehero|chegg|quizlet|studocu|reddit|youtube|wikipedia/.test(host)) score -= 10;
  return score;
}

async function searchUniversity(university) {
  const query = `"${university}" "quantum mechanics" syllabus textbook`;
  const braveUrl = `https://search.brave.com/search?q=${encodeURIComponent(query)}&source=web`;
  const braveResponse = await fetchUrl(braveUrl, {
    cacheDir: BRAVE_CACHE_DIR,
    cacheKey: `brave:${query}`,
    maxBytes: 900_000,
  });
  let parsed = parseBrave(braveResponse.body);
  if (parsed.length === 0) {
    const broaderQuery = `${university} quantum mechanics syllabus textbook`;
    const broaderUrl = `https://search.brave.com/search?q=${encodeURIComponent(broaderQuery)}&source=web`;
    const broaderResponse = await fetchUrl(broaderUrl, {
      cacheDir: BRAVE_CACHE_DIR,
      cacheKey: `brave:${broaderQuery}`,
      maxBytes: 900_000,
    });
    parsed = parseBrave(broaderResponse.body).map((result) => ({ ...result, query: broaderQuery }));
  }
  if (parsed.length > 0) {
    return parsed
      .map((result) => ({
        ...result,
        query: result.query || query,
        host: hostOf(result.url),
        score: scoreSearchResult(result),
      }))
      .sort((a, b) => b.score - a.score)
      .slice(0, 10);
  }

  const url = `https://duckduckgo.com/html/?q=${encodeURIComponent(query)}`;
  const response = await fetchUrl(url, {
    cacheDir: SEARCH_CACHE_DIR,
    cacheKey: `ddg:${query}`,
    maxBytes: 700_000,
  });
  const results = parseDdg(response.body)
    .map((result) => ({
      ...result,
      query,
      host: hostOf(result.url),
      score: scoreSearchResult(result),
    }))
    .sort((a, b) => b.score - a.score)
    .slice(0, 10);
  return results;
}

function inferSourceType(result, fetched) {
  const text = `${result.title} ${result.url} ${result.snippet}`.toLowerCase();
  const fetchedText = fetched ? cleanText(fetched.body).slice(0, 3000).toLowerCase() : "";
  if (/syllabus/.test(text) || /syllabus/.test(fetchedText)) return "official/course syllabus candidate";
  if (/reading list|readings|course reserve|leganto|talis/.test(text + " " + fetchedText)) return "reading list candidate";
  if (/opencourseware|ocw/.test(text)) return "open courseware";
  if (/course catalog|module|course outline/.test(text + " " + fetchedText)) return "course catalog/module candidate";
  return "search candidate";
}

function inferCourseTitle(result) {
  let title = result.title.replace(/\s+/g, " ").trim();
  const syllabusPipe = title.match(/Syllabus\s*\|\s*([^|]+)\s*\|/i);
  if (syllabusPipe) return syllabusPipe[1].trim();
  const beforeOcw = title.match(/^(.+?)\s*\|\s*.*OpenCourseWare/i);
  if (beforeOcw) return beforeOcw[1].trim();
  title = title.replace(/\s*-\s*MIT OpenCourseWare.*$/i, "");
  title = title.replace(/\s*\|\s*.*$/i, "");
  return title || "Unknown course";
}

function inferCourseCode(result) {
  const candidates = `${result.title} ${result.url}`;
  const mit = result.url.match(/\/courses\/([0-9]+-[0-9a-z]+)-/i);
  if (mit) return mit[1].toUpperCase();
  const code = candidates.match(/\b([A-Z]{2,5}\s?[0-9]{2,5}[A-Z]?)\b/);
  return code ? code[1].replace(/\s+/, " ") : "";
}

function inferLevel(result, courseTitle) {
  const text = `${courseTitle} ${result.title} ${result.url} ${result.snippet}`.toLowerCase();
  if (/master|msc|m\.sc|postgraduate|graduate|advanced quantum|quantum field theory|many[- ]body/.test(text)) {
    return "MSc/graduate (inferred)";
  }
  if (/bsc|b\.sc|bachelor|undergraduate|introductory|introduction|quantum physics i|quantum mechanics i/.test(text)) {
    return "BSc/undergraduate (inferred)";
  }
  return "Unknown";
}

function officialEstimate(result) {
  const host = hostOf(result.url);
  if (/\.edu$|\.edu\.|\.ac\.uk$|\.edu\.au$|\.edu\.cn$|\.ac\.jp$|\.ac\.kr$|\.edu\.sg$/.test(host)) {
    return "likely official academic domain";
  }
  if (/ocw|opencourseware|course|catalog|module|physics|university/.test(`${host} ${result.url}`.toLowerCase())) {
    return "possibly official academic page";
  }
  return "not established";
}

function pageLooksCourseLike(result, text) {
  const combined = `${result.title} ${result.url} ${result.snippet} ${text.slice(0, 3000)}`.toLowerCase();
  return /syllabus|readings|reading list|textbook|course|module|lecture/.test(combined);
}

function findBookMatches(text) {
  const found = [];
  for (const book of BOOKS) {
    const matchedPatterns = book.patterns.filter((pattern) => pattern.test(text));
    if (matchedPatterns.length > 0) {
      found.push({
        ...book,
        matched_terms: matchedPatterns.map((pattern) => pattern.source).join("; "),
      });
    }
  }
  return found;
}

function inferRequiredRecommended(text, book) {
  const idx = Math.max(
    text.toLowerCase().indexOf(book.canonical_title.toLowerCase().slice(0, 18)),
    text.toLowerCase().indexOf(book.canonical_authors.split(";")[0].toLowerCase().split(" ").pop()),
  );
  const window = idx >= 0 ? text.slice(Math.max(0, idx - 250), Math.min(text.length, idx + 250)).toLowerCase() : text.slice(0, 500).toLowerCase();
  if (/required|required text|main text|the text is|textbook for the course/.test(window)) return "required/main text (inferred)";
  if (/recommended|supplement|optional|useful|reference|further reading/.test(window)) return "recommended/supplemental (inferred)";
  if (/textbook|texts|reading|readings/.test(window)) return "textbook/reading (inferred)";
  return "matched mention";
}

function confidenceFor(result, fetched, bookMatched, courseLike) {
  if (!bookMatched) return "Candidate";
  const academic = officialEstimate(result) === "likely official academic domain";
  const fetchedOk = fetched?.status >= 200 && fetched?.status < 300;
  const content = `${result.title} ${result.url} ${result.snippet} ${fetched?.body ?? ""}`.toLowerCase();
  if (fetchedOk && academic && courseLike && /syllabus|textbook|required|readings|reading list/.test(content)) return "A";
  if (fetchedOk && academic && courseLike) return "B";
  if (academic && /textbook|isbn|required|recommended/.test(result.snippet.toLowerCase())) return "C";
  if (/textbook|isbn|required|recommended/.test(result.snippet.toLowerCase())) return "D";
  return "E";
}

async function main() {
  await fs.mkdir(DATA_DIR, { recursive: true });
  await fs.mkdir(SEARCH_CACHE_DIR, { recursive: true });
  await fs.mkdir(BRAVE_CACHE_DIR, { recursive: true });
  await fs.mkdir(PAGE_CACHE_DIR, { recursive: true });

  const ranking = await fetchArwuPayload();
  const universities = ranking.universities;
  await fs.writeFile(path.join(DATA_DIR, "universities_arwu_2025_top100.json"), JSON.stringify(universities, null, 2), "utf8");
  await writeCsv(path.join(DATA_DIR, "universities_arwu_2025_top100.csv"), universities, [
    "rank",
    "university",
    "country",
    "regional_rank",
    "ranking_source",
    "ranking_year",
    "ranking_url",
    "arwu_profile_url",
    "overall_score",
  ]);

  const candidates = [];
  const courses = [];
  const textbooks = [];
  const gaps = [];

  for (let i = 0; i < universities.length; i += 1) {
    const uni = universities[i];
    console.error(`[${i + 1}/${universities.length}] Searching ${uni.university}`);
    const results = await searchUniversity(uni.university);
    for (const [rankWithinSearch, result] of results.entries()) {
      candidates.push({
        university: uni.university,
        rank: uni.rank,
        country: uni.country,
        search_rank: rankWithinSearch + 1,
        score: result.score,
        title: result.title,
        url: result.url,
        host: result.host,
        snippet: result.snippet,
        query: result.query,
      });
    }

    const fetchable = results
      .filter((result) => result.score >= 5 && /^https?:\/\//.test(result.url))
      .slice(0, 8);

    let verifiedCount = 0;
    let candidateCount = 0;
    for (const result of fetchable) {
      await sleep(250);
      let fetched = null;
      if (!/\.pdf($|\?)/i.test(result.url)) {
        fetched = await fetchUrl(result.url, {
          cacheDir: PAGE_CACHE_DIR,
          cacheKey: `page:${result.url}`,
          maxBytes: 1_500_000,
        });
      }
      const pageText = fetched?.body ? cleanText(fetched.body) : "";
      const textForMatch = `${result.title} ${result.snippet} ${pageText}`;
      const courseTitle = inferCourseTitle(result);
      const courseCode = inferCourseCode(result);
      const level = inferLevel(result, courseTitle);
      const courseLike = pageLooksCourseLike(result, pageText);
      const sourceType = inferSourceType(result, fetched);
      const matches = findBookMatches(textForMatch);
      candidateCount += 1;

      courses.push({
        university: uni.university,
        rank: uni.rank,
        country: uni.country,
        course_code: courseCode,
        course_title: courseTitle,
        level,
        degree_context: "",
        department: "",
        academic_year: "",
        term: "",
        source_url: result.url,
        source_type: sourceType,
        source_status: matches.length ? "textbook match found" : "candidate checked; no recognized textbook match",
        confidence: matches.length ? confidenceFor(result, fetched, true, courseLike) : "Candidate",
        notes: officialEstimate(result),
      });

      for (const match of matches) {
        const confidence = confidenceFor(result, fetched, true, courseLike);
        if (["A", "B", "C"].includes(confidence)) verifiedCount += 1;
        textbooks.push({
          university: uni.university,
          rank: uni.rank,
          country: uni.country,
          course_code: courseCode,
          course_title: courseTitle,
          level,
          textbook_title: match.canonical_title,
          authors: match.canonical_authors,
          edition: "",
          publisher: "",
          year: "",
          isbn: "",
          required_or_recommended: inferRequiredRecommended(textForMatch, match),
          source_evidence_note: `Matched ${match.canonical_title}; ${sourceType}; ${officialEstimate(result)}`,
          source_url: result.url,
          confidence,
          category: match.category,
          matched_terms: match.matched_terms,
        });
      }
    }

    if (verifiedCount === 0) {
      gaps.push({
        university: uni.university,
        rank: uni.rank,
        country: uni.country,
        course_or_program: "BSc/MSc quantum courses",
        missing_item: "Publicly verified textbook assignment",
        reason:
          candidateCount === 0
            ? "No strong public search candidates in automated pass"
            : "Search candidates found, but no A-C textbook match verified by automated pass",
        attempted_sources: String(results.length),
        recommended_next_action: "Check official department catalogue, reading-list system, bookstore by course code, then contact department/instructor.",
      });
    }
    await sleep(700);
  }

  const canonicalMap = new Map();
  for (const row of textbooks) {
    const key = `${row.textbook_title}--${row.authors}`;
    const current = canonicalMap.get(key) ?? {
      canonical_title: row.textbook_title,
      canonical_authors: row.authors,
      category: row.category,
      used_by_university_count: 0,
      used_in_course_count: 0,
      bsc_count: 0,
      msc_count: 0,
      unknown_level_count: 0,
      confidence_a_to_c_count: 0,
      source_urls: new Set(),
      universities: new Set(),
      courses: new Set(),
    };
    current.universities.add(row.university);
    current.courses.add(`${row.university}::${row.course_code}::${row.course_title}`);
    current.source_urls.add(row.source_url);
    if (/BSc|undergraduate/.test(row.level)) current.bsc_count += 1;
    else if (/MSc|graduate/.test(row.level)) current.msc_count += 1;
    else current.unknown_level_count += 1;
    if (["A", "B", "C"].includes(row.confidence)) current.confidence_a_to_c_count += 1;
    canonicalMap.set(key, current);
  }
  const canonicalBooks = Array.from(canonicalMap.values())
    .map((row) => ({
      canonical_title: row.canonical_title,
      canonical_authors: row.canonical_authors,
      category: row.category,
      used_by_university_count: row.universities.size,
      used_in_course_count: row.courses.size,
      bsc_count: row.bsc_count,
      msc_count: row.msc_count,
      unknown_level_count: row.unknown_level_count,
      confidence_a_to_c_count: row.confidence_a_to_c_count,
      source_url_count: row.source_urls.size,
    }))
    .sort((a, b) => b.used_by_university_count - a.used_by_university_count || a.canonical_title.localeCompare(b.canonical_title));

  const summary = {
    generated_at: new Date().toISOString(),
    ranking_source: "ARWU 2025",
    ranking_url: RANKING_PAGE,
    arwu_payload_url: ranking.payloadUrl,
    university_count: universities.length,
    candidate_count: candidates.length,
    course_candidate_count: courses.length,
    textbook_match_count: textbooks.length,
    verified_textbook_match_count: textbooks.filter((row) => ["A", "B", "C"].includes(row.confidence)).length,
    universities_with_verified_matches: new Set(textbooks.filter((row) => ["A", "B", "C"].includes(row.confidence)).map((row) => row.university)).size,
    universities_with_gaps: gaps.length,
  };

  await fs.writeFile(path.join(DATA_DIR, "search_candidates.json"), JSON.stringify(candidates, null, 2), "utf8");
  await fs.writeFile(path.join(DATA_DIR, "courses.json"), JSON.stringify(courses, null, 2), "utf8");
  await fs.writeFile(path.join(DATA_DIR, "textbooks.json"), JSON.stringify(textbooks, null, 2), "utf8");
  await fs.writeFile(path.join(DATA_DIR, "canonical_books.json"), JSON.stringify(canonicalBooks, null, 2), "utf8");
  await fs.writeFile(path.join(DATA_DIR, "gaps.json"), JSON.stringify(gaps, null, 2), "utf8");
  await fs.writeFile(path.join(DATA_DIR, "summary.json"), JSON.stringify(summary, null, 2), "utf8");

  await writeCsv(path.join(DATA_DIR, "search_candidates.csv"), candidates, [
    "university",
    "rank",
    "country",
    "search_rank",
    "score",
    "title",
    "url",
    "host",
    "snippet",
    "query",
  ]);
  await writeCsv(path.join(DATA_DIR, "courses.csv"), courses, [
    "university",
    "rank",
    "country",
    "course_code",
    "course_title",
    "level",
    "degree_context",
    "department",
    "academic_year",
    "term",
    "source_url",
    "source_type",
    "source_status",
    "confidence",
    "notes",
  ]);
  await writeCsv(path.join(DATA_DIR, "textbooks.csv"), textbooks, [
    "university",
    "rank",
    "country",
    "course_code",
    "course_title",
    "level",
    "textbook_title",
    "authors",
    "edition",
    "publisher",
    "year",
    "isbn",
    "required_or_recommended",
    "source_evidence_note",
    "source_url",
    "confidence",
    "category",
    "matched_terms",
  ]);
  await writeCsv(path.join(DATA_DIR, "canonical_books.csv"), canonicalBooks, [
    "canonical_title",
    "canonical_authors",
    "category",
    "used_by_university_count",
    "used_in_course_count",
    "bsc_count",
    "msc_count",
    "unknown_level_count",
    "confidence_a_to_c_count",
    "source_url_count",
  ]);
  await writeCsv(path.join(DATA_DIR, "gaps.csv"), gaps, [
    "university",
    "rank",
    "country",
    "course_or_program",
    "missing_item",
    "reason",
    "attempted_sources",
    "recommended_next_action",
  ]);

  console.error(JSON.stringify(summary, null, 2));
}

main().catch((error) => {
  console.error(error.stack || error.message);
  process.exit(1);
});
