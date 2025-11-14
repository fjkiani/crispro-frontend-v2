# ⚔️ JR1 PRE-SEEDING CHECKLIST ⚔️

**Owner:** Agent JR1  
**Mission Commander:** Zo  
**Status:** 🔄 **PRE-FLIGHT CHECK BEFORE SEEDING**  
**Goal:** Ensure database schema + parsing logic ready for 200 trial ingestion

---

## 🎯 **CRITICAL QUESTIONS BEFORE SEEDING**

### **1. AstraDB SCHEMA VERIFICATION**

**Question:** Does the AstraDB `trials_2024` table have ALL fields we need for JR2's GTM?

**Required Fields (for clinical use - Ayesha):**
- ✅ `trial_id` (NCT ID)
- ✅ `title`
- ✅ `phase`
- ✅ `status` (Recruiting/Not yet recruiting)
- ✅ `interventions` (drug names, classes)
- ✅ `inclusion_criteria` (text)
- ✅ `exclusion_criteria` (text)
- ✅ `locations` (sites)
- ✅ `disease` (ovarian cancer)
- ✅ `treatment_line` (frontline)

**NEW FIELDS NEEDED (for GTM use - JR2):**
- ❓ `sponsor_name` (e.g., "AstraZeneca", "Genentech")
- ❓ `sponsor_contact_email` (if available from ClinicalTrials.gov)
- ❓ `principal_investigator_name` (lead PI)
- ❓ `pi_contact_email` (if available)
- ❓ `study_coordinator_email` (if available)
- ❓ `primary_endpoint` (e.g., "PFS", "ORR")
- ❓ `mechanism_tags` (e.g., ["PARPi", "anti-angiogenic", "immunotherapy"])
- ❓ `biomarker_requirements` (e.g., "HRD", "BRCA", "PDL1")
- ❓ `site_count` (number of locations)
- ❓ `estimated_enrollment` (target n)

**ACTION FOR JR1:**
1. ✅ Review existing `trials_2024` schema in AstraDB
2. ❓ If missing GTM fields, add them to schema BEFORE seeding
3. ❓ Update parsing logic in `api/services/hybrid_trial_search.py` or wherever trials are ingested

---

### **2. PARSING LOGIC FOR GTM FIELDS**

**Question:** Where do we get sponsor/PI/contact info from?

**Primary Source:** ClinicalTrials.gov API (`https://clinicaltrials.gov/api/v2/studies`)

**Key API Endpoints:**
- `/studies/{NCT_ID}` → Returns full trial JSON with:
  - `protocolSection.sponsorCollaboratorsModule.leadSponsor.name`
  - `protocolSection.sponsorCollaboratorsModule.leadSponsor.class` (INDUSTRY, NIH, OTHER)
  - `protocolSection.contactsLocationsModule.centralContacts` (study-wide contact)
  - `protocolSection.contactsLocationsModule.locations[*].contacts` (site-level contacts)
  - `protocolSection.outcomesModule.primaryOutcomes[0].measure` (primary endpoint)
  - `protocolSection.designModule.enrollmentInfo.count` (target enrollment)

**Example Parsing Code (for JR1 to add):**
```python
def parse_gtm_fields(trial_json):
    """Extract GTM-specific fields from ClinicalTrials.gov API response"""
    protocol = trial_json.get('protocolSection', {})
    
    # Sponsor
    sponsor_module = protocol.get('sponsorCollaboratorsModule', {})
    sponsor_name = sponsor_module.get('leadSponsor', {}).get('name', 'Unknown')
    
    # Contacts
    contacts_module = protocol.get('contactsLocationsModule', {})
    central_contacts = contacts_module.get('centralContacts', [])
    sponsor_contact_email = central_contacts[0].get('email', None) if central_contacts else None
    
    # PI (from first location)
    locations = contacts_module.get('locations', [])
    pi_name = None
    pi_email = None
    if locations:
        first_location_contacts = locations[0].get('contacts', [])
        if first_location_contacts:
            pi_name = first_location_contacts[0].get('name', None)
            pi_email = first_location_contacts[0].get('email', None)
    
    # Primary endpoint
    outcomes_module = protocol.get('outcomesModule', {})
    primary_outcomes = outcomes_module.get('primaryOutcomes', [])
    primary_endpoint = primary_outcomes[0].get('measure', 'Not specified') if primary_outcomes else 'Not specified'
    
    # Enrollment
    design_module = protocol.get('designModule', {})
    enrollment_info = design_module.get('enrollmentInfo', {})
    estimated_enrollment = enrollment_info.get('count', 0)
    
    # Site count
    site_count = len(locations)
    
    return {
        'sponsor_name': sponsor_name,
        'sponsor_contact_email': sponsor_contact_email,
        'principal_investigator_name': pi_name,
        'pi_contact_email': pi_email,
        'primary_endpoint': primary_endpoint,
        'estimated_enrollment': estimated_enrollment,
        'site_count': site_count
    }
```

**ACTION FOR JR1:**
- ✅ Add this parsing logic to trial ingestion script
- ✅ Test on 1-2 trials before mass seeding
- ✅ Handle missing fields gracefully (None/null, not empty strings)

---

### **3. MECHANISM TAGGING**

**Question:** How do we auto-tag trials with mechanism types for JR2's 1-pagers?

**Approach:** Simple keyword matching on `interventions` + `title`

**Example Logic:**
```python
def tag_mechanisms(trial_title, interventions):
    """Auto-tag trial mechanisms based on interventions"""
    mechanisms = []
    
    # Combine title + intervention names
    text = f"{trial_title} {' '.join(interventions)}".lower()
    
    # PARPi
    if any(word in text for word in ['parp', 'olaparib', 'niraparib', 'rucaparib', 'talazoparib']):
        mechanisms.append('PARPi')
    
    # Anti-angiogenic
    if any(word in text for word in ['bevacizumab', 'avastin', 'anti-vegf', 'angiogenesis']):
        mechanisms.append('anti-angiogenic')
    
    # Immunotherapy
    if any(word in text for word in ['pembrolizumab', 'nivolumab', 'checkpoint', 'pd-1', 'pd-l1', 'immunotherapy']):
        mechanisms.append('immunotherapy')
    
    # Chemotherapy
    if any(word in text for word in ['carboplatin', 'paclitaxel', 'docetaxel', 'gemcitabine', 'platinum']):
        mechanisms.append('chemotherapy')
    
    # ADC
    if any(word in text for word in ['antibody-drug conjugate', 'adc', 'mirvetuximab', 'sacituzumab']):
        mechanisms.append('ADC')
    
    return mechanisms if mechanisms else ['other']
```

**ACTION FOR JR1:**
- ✅ Add this tagging logic
- ✅ Store as JSON array in AstraDB: `mechanism_tags: ["PARPi", "anti-angiogenic"]`

---

### **4. BIOMARKER REQUIREMENTS EXTRACTION**

**Question:** How do we extract biomarker requirements for JR2's 1-pagers?

**Approach:** Parse eligibility criteria for common biomarkers

**Example Logic:**
```python
def extract_biomarker_requirements(inclusion_criteria):
    """Extract required biomarkers from eligibility text"""
    biomarkers = []
    text = inclusion_criteria.lower()
    
    # HRD
    if any(word in text for word in ['hrd', 'homologous recombination deficiency', 'genomic instability']):
        biomarkers.append('HRD')
    
    # BRCA
    if any(word in text for word in ['brca1', 'brca2', 'brca mutation', 'germline brca']):
        biomarkers.append('BRCA')
    
    # PDL1
    if any(word in text for word in ['pd-l1', 'pdl1', 'pd-l1 expression']):
        biomarkers.append('PDL1')
    
    # TMB
    if any(word in text for word in ['tmb', 'tumor mutational burden']):
        biomarkers.append('TMB')
    
    # MSI
    if any(word in text for word in ['msi-h', 'microsatellite instability', 'dmmr']):
        biomarkers.append('MSI-H')
    
    return biomarkers if biomarkers else None
```

**ACTION FOR JR1:**
- ✅ Add this extraction logic
- ✅ Store as JSON array: `biomarker_requirements: ["HRD", "BRCA"]` or `null`

---

## 📊 **UPDATED TRIAL SCHEMA (FINAL)**

**AstraDB `trials_2024` Table Schema:**

```json
{
    "trial_id": "NCT05551052",
    "title": "Study of Niraparib + Bevacizumab in First-Line Ovarian Cancer",
    "phase": "Phase 3",
    "status": "Recruiting",
    "interventions": ["Niraparib", "Bevacizumab"],
    "inclusion_criteria": "...",
    "exclusion_criteria": "...",
    "locations": [...],
    "disease": "ovarian_cancer",
    "treatment_line": "first-line",
    
    // NEW GTM FIELDS
    "sponsor_name": "GlaxoSmithKline",
    "sponsor_contact_email": "clinical.trials@gsk.com",
    "principal_investigator_name": "Dr. Jane Smith",
    "pi_contact_email": "jsmith@msk.org",
    "study_coordinator_email": null,
    "primary_endpoint": "Progression-Free Survival (PFS)",
    "mechanism_tags": ["PARPi", "anti-angiogenic"],
    "biomarker_requirements": ["HRD"],
    "site_count": 47,
    "estimated_enrollment": 538
}
```

---

## ✅ **PRE-SEEDING CHECKLIST FOR JR1**

**Before running mass seeding:**

- [ ] **Schema Update**: Add GTM fields to AstraDB `trials_2024` table
- [ ] **Parsing Logic**: Implement `parse_gtm_fields()` function
- [ ] **Mechanism Tagging**: Implement `tag_mechanisms()` function
- [ ] **Biomarker Extraction**: Implement `extract_biomarker_requirements()` function
- [ ] **Test Run**: Seed 5-10 trials manually, verify all fields populate correctly
- [ ] **Error Handling**: Graceful fallbacks for missing contacts (use `null`, not empty strings)
- [ ] **Coordinate with Zo**: Confirm schema changes won't break existing Ayesha trial search

**Once all checkboxes complete:**
✅ **READY TO SEED 200 TRIALS** → Proceed with mass ingestion

---

## 🎯 **SUCCESS CRITERIA**

**After seeding, verify:**
- ✅ All 200 trials have `sponsor_name` (0% null)
- ✅ At least 60% have `sponsor_contact_email` or `pi_contact_email`
- ✅ At least 80% have `mechanism_tags` (not "other")
- ✅ At least 50% have `biomarker_requirements` (non-null)
- ✅ All trials have `primary_endpoint`, `site_count`, `estimated_enrollment`

**If success criteria met:**
→ **JR2 can begin GTM automation** (mass 1-pager generation)

**If success criteria NOT met:**
→ **JR1 must re-parse** with improved extraction logic

---

**STATUS:** 🔄 **AWAITING JR1 EXECUTION**

