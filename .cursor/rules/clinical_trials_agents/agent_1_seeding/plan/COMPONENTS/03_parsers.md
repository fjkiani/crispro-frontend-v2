# 📝 MODULE 3: PARSERS

## **Purpose**
Parse ClinicalTrials.gov API v2 responses into our schema format

## **File Locations**
- `scripts/agent_1_seeding/parsers/study_parser.py` - Main orchestrator
- `scripts/agent_1_seeding/parsers/biomarker_extractor.py` - Biomarker extraction
- `scripts/agent_1_seeding/parsers/locations_parser.py` - Locations extraction

## **Key Functions**

### **`parse_ctgov_study(study: Dict) -> Dict`** (study_parser.py)
Main parser that orchestrates sub-parsers:
- Extracts IDs, status, phase, description
- Parses inclusion/exclusion criteria
- Calls `extract_biomarkers()` and `parse_locations_data()`
- Returns complete schema dict with new fields

### **`extract_biomarkers(eligibility_text: str) -> List[str]`** (biomarker_extractor.py)
Keyword-based biomarker extraction:
- Uses `BIOMARKER_KEYWORDS` from config
- Case-insensitive matching
- Returns deduplicated list

### **`parse_locations_data(study: Dict) -> List[Dict]`** (locations_parser.py)
Locations data extraction with error handling:
- Extracts facility, city, state, zip, contacts
- Validates critical fields (facility/state)
- Skips locations missing critical fields
- Returns JSON-serializable list

## **Error Handling**
- ✅ Check `len(contacts) > 0` before indexing
- ✅ Skip locations missing facility/state
- ✅ Graceful handling of missing API fields

## **Dependencies**
- `config.py` - Biomarker keywords, disease categories

## **Acceptance Criteria**
- [ ] Parser handles all API v2 fields
- [ ] Biomarker extraction works (BRCA1/2, HRD)
- [ ] Locations data properly formatted as JSON
- [ ] Error handling for missing/empty data
- [ ] KM/HRD updates (cloned/copied field names) reflected

## **Time Estimate:** 1 hour

