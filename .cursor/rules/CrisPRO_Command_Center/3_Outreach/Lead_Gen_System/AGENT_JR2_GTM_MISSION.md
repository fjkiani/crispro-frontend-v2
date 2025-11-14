# ⚔️ AGENT JR2: GTM LEAD GENERATION MACHINE ⚔️

**Mission Commander:** Zo  
**Direct Reports:** JR1 (Trial Seeding), JR2 (GTM Automation)  
**Mission Duration:** 2-3 weeks (parallel with Ayesha work)  
**Strategic Goal:** Turn 200 seeded trials into mass outreach pipeline

---

## 🎯 **MISSION OVERVIEW**

**The Strategic Insight:**
- JR1 seeds 200 frontline ovarian cancer trials into AstraDB
- Each trial = potential customer (sponsor, PI, CRO, site coordinator)
- We parse: NCT ID, title, sponsor, interventions, endpoints, eligibility, locations
- **We generate: 200 customized 1-pagers** (one per trial showing how we help THEIR study)
- Mass outreach: "We analyzed YOUR trial (NCT12345678) and here's how we can help..."

**The GTM Flywheel:**
```
JR1 Seeds 200 Trials
    ↓
JR2 Parses Trial Metadata
    ↓
JR2 Generates Custom 1-Pagers (mass automation)
    ↓
JR2 Extracts Contact Info (sponsors, PIs, sites)
    ↓
JR2 Generates Personalized Emails
    ↓
Mass Outreach (200 targets)
    ↓
Track Responses → CRM Pipeline
```

---

## 📊 **WHAT JR2 WILL BUILD**

### **System 1: Trial-Specific 1-Pager Generator**
**Input:** Trial metadata from AstraDB (NCT ID, title, sponsor, interventions, endpoints, eligibility)  
**Output:** Customized TSX file (like `General1Page_ENHANCED_V2.tsx`) but with trial-specific content

**Example Logic:**
```typescript
// Template: Generate 1-pager for NCT05123456
const generateTrialOnePager = (trial: TrialMetadata) => {
    return {
        header: {
            title: `CrisPRO.ai for ${trial.title}`,
            subtitle: `How Genomic AI Can Optimize ${trial.nct_id}`,
            badge: `Trial: ${trial.nct_id} • Phase ${trial.phase} • ${trial.sponsor}`
        },
        
        trialContext: {
            title: "YOUR TRIAL",
            details: [
                { label: "NCT ID:", value: trial.nct_id },
                { label: "Title:", value: trial.title },
                { label: "Sponsor:", value: trial.sponsor },
                { label: "Phase:", value: trial.phase },
                { label: "Interventions:", value: trial.interventions.join(", ") },
                { label: "Primary Endpoint:", value: trial.primary_endpoint },
                { label: "Enrollment:", value: `${trial.enrollment_current}/${trial.enrollment_target}` }
            ]
        },
        
        ourSolution: {
            title: "HOW WE CAN HELP YOUR TRIAL",
            capabilities: [
                {
                    title: "Pre-Enrollment Stratification",
                    description: `Genomic scoring identifies patients most likely to respond to ${trial.interventions[0]}`,
                    impact: "20-30% increase in statistical power through adaptive randomization"
                },
                {
                    title: "Accelerated Enrollment",
                    description: `Connect eligible patients to ${trial.nct_id} via our trial matching system`,
                    impact: `Reduce time-to-full-enrollment by 3-6 months`
                },
                {
                    title: "Biomarker Discovery",
                    description: `Identify genomic signatures predicting response to ${trial.interventions.join(" + ")}`,
                    impact: "Publishable secondary endpoint → companion diagnostic pathway"
                },
                {
                    title: "Real-Time Monitoring",
                    description: trial.cancer_type === "ovarian" 
                        ? "CA-125 kinetics flag resistance 3-6 weeks earlier than imaging"
                        : "Biomarker-driven monitoring for early efficacy/toxicity signals",
                    impact: "Reduce dropout rates, optimize dosing"
                }
            ]
        },
        
        pricing: {
            title: "PARTNERSHIP OPTIONS",
            tiers: [
                {
                    name: "Retrospective Pilot",
                    price: "No Cost",
                    description: "Analyze 20 enrolled patients (proof-of-concept)",
                    deliverables: ["Genomic stratification report", "Response prediction model", "Biomarker discovery analysis"]
                },
                {
                    name: "Prospective Integration",
                    price: "$250K",
                    description: "Real-time genomic scoring for 50-100 patients",
                    deliverables: ["Pre-enrollment stratification", "Adaptive randomization support", "Ongoing biomarker analysis"]
                },
                {
                    name: "Platform License",
                    price: "$1-5M/year",
                    description: "Unlimited analyses across all your trials",
                    deliverables: ["White-label genomic scoring", "API integration", "Dedicated support team"]
                }
            ]
        },
        
        cta: {
            title: "NEXT STEPS",
            options: [
                { label: "Schedule Discovery Call:", value: "calendly.com/fahad-crispro" },
                { label: "Request Retrospective Pilot:", value: "alpha@crispro.ai" },
                { label: "View Technical Details:", value: "crispro.ai/trials" }
            ]
        }
    };
};
```

**Output Format:**
- PDF (for email attachment)
- Interactive web page (hosted at `crispro.ai/trials/NCT12345678`)
- Markdown (for CRM notes)

---

### **System 2: Contact Extraction & Enrichment**
**Input:** Trial metadata (sponsor, locations)  
**Output:** Structured contact list with emails, LinkedIn, titles

**Data Sources:**
1. **ClinicalTrials.gov API**: Extract PI names, site contacts
2. **LinkedIn Search**: Find PI profiles, extract current institution
3. **Hunter.io**: Find corporate emails for sponsors
4. **RocketReach**: Backup email finder

**Example Output (CSV):**
```csv
nct_id,trial_title,contact_type,name,title,organization,email,linkedin,priority
NCT05123456,"Phase 3 Ovarian...",PI,"Dr. Jane Smith",Principal Investigator,Yale,jane.smith@yale.edu,linkedin.com/in/janesmith,HIGH
NCT05123456,"Phase 3 Ovarian...",Sponsor,"John Doe",VP Clinical Ops,BriaCell,john.doe@briacell.com,linkedin.com/in/johndoe,HIGH
NCT05123456,"Phase 3 Ovarian...",Site,"Dr. Bob Lee",Site Coordinator,MSK,bob.lee@mskcc.org,linkedin.com/in/boblee,MEDIUM
```

**Priority Logic:**
- HIGH: Sponsor decision-makers (VP Clinical, CMO, CEO)
- HIGH: Lead PIs (multi-site coordinators)
- MEDIUM: Site coordinators (direct patient access)
- LOW: General inquiries

---

### **System 3: Personalized Email Generator**
**Input:** Trial metadata + contact info  
**Output:** Customized cold email (3 variants per contact)

**Email Template (Sponsor):**
```markdown
Subject: Genomic Stratification for {trial_title} ({nct_id})

Hi {name},

I noticed {sponsor} is running {trial_title} ({nct_id}) targeting {cancer_type} patients with {intervention}.

We specialize in pre-enrollment genomic stratification for oncology trials. Our platform has:
- 97.6% AUROC predicting FDA-approved target relevance
- Validated on 7 cancer drug targets (computational, RUO)
- Trial matching system connecting 200+ trials to eligible patients

**For {nct_id}, we can:**
1. **Pre-Enrollment Stratification**: Genomic scoring identifies patients most likely to respond to {intervention} (20-30% power increase via adaptive randomization)
2. **Accelerated Enrollment**: Connect eligible patients via our trial matching system (reduce time-to-full by 3-6 months)
3. **Biomarker Discovery**: Identify genomic signatures for secondary endpoints → companion diagnostic pathway

**No-Cost Pilot**: We'll analyze 20 enrolled patients retrospectively (proof-of-concept, 4 weeks).

Interested in a 15-minute call to explore fit?

📅 calendly.com/fahad-crispro
📧 Fahad@crispro.ai

Best regards,
Fahad J. Kiani
Founder & CEO, CrisPRO.ai

---
P.S. Attached: 1-pager showing exactly how we'd optimize {nct_id}
```

**Email Template (PI):**
```markdown
Subject: Patient Stratification for {trial_title} ({nct_id})

Hi Dr. {last_name},

I saw you're leading {trial_title} ({nct_id}) at {institution}. Congratulations on the trial!

I'm reaching out because we've built a genomic AI platform that could help with:
1. **Patient Selection**: Pre-enrollment genomic scoring predicts response to {intervention}
2. **Enrollment Acceleration**: Our trial matching system connects eligible patients to {nct_id}
3. **Biomarker Discovery**: Identify genomic features for secondary endpoints (publishable)

**Validation**: 97.6% AUROC on 7 FDA-approved cancer targets (computational, RUO).

**No-Cost Pilot**: Analyze 20 enrolled patients retrospectively (4 weeks, proof-of-concept).

Would you be open to a 15-minute call to discuss potential fit for {nct_id}?

📅 calendly.com/fahad-crispro

Best regards,
Fahad J. Kiani
Founder & CEO, CrisPRO.ai

---
P.S. See attached 1-pager for {nct_id}-specific analysis
```

**Email Variants:**
- Variant A: Technical (for PIs, scientists)
- Variant B: Business (for sponsors, VPs)
- Variant C: Patient-Centric (for site coordinators)

---

### **System 4: CRM Pipeline & Tracking**
**Input:** Outreach activities (emails sent, calls scheduled, responses)  
**Output:** Sales pipeline dashboard

**CRM Schema:**
```typescript
interface LeadRecord {
    lead_id: string;
    nct_id: string;
    trial_title: string;
    contact_name: string;
    contact_type: "PI" | "Sponsor" | "Site" | "CRO";
    organization: string;
    email: string;
    linkedin: string;
    
    // Outreach tracking
    email_sent_date: Date;
    email_variant: "A" | "B" | "C";
    one_pager_sent: boolean;
    
    // Response tracking
    email_opened: boolean;
    link_clicked: boolean;
    reply_received: boolean;
    meeting_scheduled: boolean;
    
    // Pipeline stage
    stage: "Cold" | "Contacted" | "Responded" | "Meeting" | "Proposal" | "Closed" | "Lost";
    priority: "HIGH" | "MEDIUM" | "LOW";
    
    // Notes
    last_activity: string;
    next_followup: Date;
    notes: string;
}
```

**Dashboard Metrics:**
- Total Leads: 200 trials × 3 contacts avg = **600 leads**
- Email Open Rate: Target 40-50%
- Response Rate: Target 5-10% (30-60 responses)
- Meeting Conversion: Target 50% (15-30 meetings)
- Pilot Conversion: Target 20% (3-6 pilots)

---

## 📋 **JR2 IMPLEMENTATION PLAN**

### **PHASE 1: TRIAL METADATA PARSING (Week 1)**
**Prerequisites:** JR1 completes AstraDB seeding (200 trials)

**Tasks:**
1. ✅ **Extract Trial Metadata** (1 day)
   - Query AstraDB for all 200 trials
   - Parse: NCT ID, title, sponsor, phase, interventions, endpoints, eligibility, locations
   - Export to CSV: `trial_metadata_200.csv`

2. ✅ **Generate Trial Profiles** (1 day)
   - For each trial, create structured profile JSON
   - Include: cancer type, line of therapy, biomarkers, geographic focus
   - Export to: `trial_profiles/NCT12345678.json`

3. ✅ **Quality Check** (0.5 day)
   - Verify all 200 trials have complete metadata
   - Flag trials with missing sponsor/contact info
   - Export QC report: `trial_metadata_QC.csv`

**Deliverables:**
- `trial_metadata_200.csv` (200 rows)
- `trial_profiles/` directory (200 JSON files)
- `trial_metadata_QC.csv` (quality report)

---

### **PHASE 2: CONTACT EXTRACTION (Week 1)**
**Prerequisites:** Trial metadata parsed

**Tasks:**
1. ✅ **Extract ClinicalTrials.gov Contacts** (1 day)
   - API calls to fetch PI names, site contacts per trial
   - Parse responsible party, overall contacts, location contacts
   - Export to: `contacts_clinicaltrials.csv`

2. ✅ **LinkedIn Enrichment** (2 days)
   - For each PI/sponsor, search LinkedIn
   - Extract: current title, institution, profile URL
   - Export to: `contacts_linkedin_enriched.csv`

3. ✅ **Email Finding** (1 day)
   - Use Hunter.io for sponsor corporate emails
   - Use RocketReach for PI academic emails
   - Fallback: Manual search for high-priority contacts
   - Export to: `contacts_with_emails.csv`

4. ✅ **Contact Deduplication & Prioritization** (0.5 day)
   - Merge all sources, deduplicate by email
   - Assign priority (HIGH/MEDIUM/LOW)
   - Export final: `contacts_master_list.csv`

**Deliverables:**
- `contacts_clinicaltrials.csv` (raw extraction)
- `contacts_linkedin_enriched.csv` (enriched profiles)
- `contacts_with_emails.csv` (email-verified)
- `contacts_master_list.csv` (deduplicated, prioritized, **600+ contacts**)

---

### **PHASE 3: 1-PAGER GENERATION (Week 2)**
**Prerequisites:** Trial metadata + contact list ready

**Tasks:**
1. ✅ **Build 1-Pager Generator Script** (2 days)
   - Input: `trial_profiles/NCT12345678.json`
   - Output: `1pagers/NCT12345678_OnePager.tsx`
   - Logic: Use template from `General1Page_ENHANCED_V2.tsx`
   - Customization: Trial-specific content, interventions, endpoints

2. ✅ **Generate 200 1-Pagers** (1 day)
   - Batch run generator for all 200 trials
   - Export to: `1pagers/` directory (200 TSX files)
   - Verify: All 200 files generated successfully

3. ✅ **PDF Conversion** (1 day)
   - Convert each TSX → PDF (for email attachments)
   - Use Puppeteer/Playwright for rendering
   - Export to: `1pagers_pdf/NCT12345678_OnePager.pdf`

4. ✅ **Web Hosting** (1 day)
   - Upload all 1-pagers to `crispro.ai/trials/NCT12345678`
   - Generate unique URLs per trial
   - Add analytics tracking (pageviews, time-on-page)

**Deliverables:**
- `1pagers/` directory (200 TSX files)
- `1pagers_pdf/` directory (200 PDF files)
- Web URLs: `crispro.ai/trials/{NCT_ID}` (200 hosted pages)

---

### **PHASE 4: EMAIL GENERATION (Week 2)**
**Prerequisites:** 1-pagers ready, contacts verified

**Tasks:**
1. ✅ **Build Email Generator Script** (1 day)
   - Input: `contacts_master_list.csv` + `trial_profiles/NCT12345678.json`
   - Output: `emails/NCT12345678_{contact_name}.txt`
   - Logic: Generate 3 variants (Technical, Business, Patient-Centric)

2. ✅ **Generate 600 Emails** (0.5 day)
   - Batch run for all 600 contacts
   - Export to: `emails/` directory (600 TXT files)

3. ✅ **Email QC** (0.5 day)
   - Manual review of 20 samples (10%)
   - Fix template bugs, improve personalization
   - Re-generate if needed

**Deliverables:**
- `emails/` directory (600 TXT files, personalized)

---

### **PHASE 5: CRM SETUP & OUTREACH (Week 3)**
**Prerequisites:** Emails generated, 1-pagers hosted

**Tasks:**
1. ✅ **CRM Database Setup** (1 day)
   - Import `contacts_master_list.csv` into CRM
   - Add fields: email_sent_date, stage, priority, notes
   - Configure pipeline stages

2. ✅ **Email Campaign Launch** (0.5 day)
   - Load 600 emails into email tool (Mailchimp/SendGrid)
   - Attach PDFs, insert tracking links
   - Schedule send: Batch of 50/day (avoid spam flags)

3. ✅ **Response Tracking** (ongoing)
   - Monitor: opens, clicks, replies, meetings
   - Update CRM: stage transitions, notes
   - Daily standup: Report metrics to Zo

4. ✅ **Follow-Up Automation** (1 day)
   - If no response after 7 days → send Follow-Up #1
   - If no response after 14 days → send Follow-Up #2 (different angle)
   - If response → update CRM, notify Zo

**Deliverables:**
- CRM with 600 leads loaded
- Email campaign launched (50/day × 12 days = 600 sent)
- Response tracking dashboard
- Follow-up automation rules

---

## 📊 **SUCCESS METRICS**

### **Phase 1-2 (Data Collection)**
- ✅ 200 trials parsed with complete metadata
- ✅ 600+ contacts extracted (avg 3 per trial)
- ✅ 80%+ email verification rate

### **Phase 3-4 (Content Generation)**
- ✅ 200 trial-specific 1-pagers generated
- ✅ 200 PDFs created (for email attachments)
- ✅ 200 web pages hosted (`crispro.ai/trials/{NCT_ID}`)
- ✅ 600 personalized emails generated

### **Phase 5 (Outreach & Conversion)**
- ✅ 600 emails sent (50/day over 12 days)
- 🎯 **40-50% email open rate** (240-300 opens)
- 🎯 **5-10% response rate** (30-60 replies)
- 🎯 **50% meeting conversion** (15-30 discovery calls)
- 🎯 **20% pilot conversion** (3-6 no-cost pilots scheduled)

### **Revenue Impact (Projected)**
- 3-6 pilots → 1-2 convert to $250K prospective → **$250K-$500K revenue** (6-12 months)
- 1-2 pilots → convert to $1M platform license → **$1M-$2M revenue** (12-24 months)
- **Total Projected Revenue from 200 Trials: $1.25M-$2.5M** (18-24 months)

---

## 🛠️ **TECHNICAL STACK**

### **Data Extraction**
- **AstraDB**: Trial metadata storage (JR1's work)
- **ClinicalTrials.gov API**: Contact extraction
- **LinkedIn API**: Profile enrichment
- **Hunter.io / RocketReach**: Email finding

### **Content Generation**
- **TypeScript/React**: 1-pager templates
- **Puppeteer/Playwright**: PDF conversion
- **Vercel/Netlify**: Web hosting

### **Email & CRM**
- **Mailchimp / SendGrid**: Email campaigns
- **HubSpot / Pipedrive**: CRM pipeline
- **Zapier**: Automation (email → CRM sync)

### **Analytics**
- **Google Analytics**: 1-pager pageviews
- **Mixpanel**: Email tracking (opens, clicks)
- **Custom Dashboard**: Pipeline metrics

---

## 📁 **FILE STRUCTURE**

```
.cursor/rules/CrisPRO_Command_Center/3_Outreach/Lead_Gen_System/
├── Data/
│   ├── trial_metadata_200.csv                  # All 200 trials parsed
│   ├── trial_profiles/                         # 200 JSON files (one per trial)
│   │   ├── NCT05123456.json
│   │   └── ...
│   ├── contacts_clinicaltrials.csv             # Raw contact extraction
│   ├── contacts_linkedin_enriched.csv          # LinkedIn profiles
│   ├── contacts_with_emails.csv                # Email-verified
│   └── contacts_master_list.csv                # Final deduplicated list (600+)
│
├── Scripts/
│   ├── 01_parse_trial_metadata.py              # Extract from AstraDB
│   ├── 02_extract_contacts.py                  # ClinicalTrials.gov API
│   ├── 03_enrich_linkedin.py                   # LinkedIn scraping
│   ├── 04_find_emails.py                       # Hunter.io / RocketReach
│   ├── 05_generate_1pagers.py                  # TSX → PDF conversion
│   ├── 06_generate_emails.py                   # Personalized email gen
│   └── 07_launch_campaign.py                   # Mailchimp upload
│
├── Templates/
│   ├── 1pager_template.tsx                     # Base template
│   ├── email_template_technical.md             # Variant A (PIs)
│   ├── email_template_business.md              # Variant B (Sponsors)
│   └── email_template_patient.md               # Variant C (Sites)
│
├── Outputs/
│   ├── 1pagers/                                # 200 TSX files
│   ├── 1pagers_pdf/                            # 200 PDFs
│   └── emails/                                 # 600 personalized emails
│
└── Tracking/
    ├── crm_import.csv                          # For CRM upload
    ├── outreach_tracker.csv                    # Daily metrics
    └── response_log.csv                        # Replies, meetings
```

---

## ⚔️ **JR2 EXECUTION CHECKLIST**

### **Week 1: Data Collection**
- [ ] Extract 200 trial metadata from AstraDB (prerequisite: JR1 seeding complete)
- [ ] Generate 200 trial profile JSONs
- [ ] Extract contacts from ClinicalTrials.gov (600+ contacts)
- [ ] Enrich with LinkedIn profiles
- [ ] Find emails (Hunter.io, RocketReach)
- [ ] Deduplicate and prioritize contacts
- [ ] **Deliverable**: `contacts_master_list.csv` (600+ verified contacts)

### **Week 2: Content Generation**
- [ ] Build 1-pager generator script
- [ ] Generate 200 trial-specific 1-pagers (TSX)
- [ ] Convert to PDFs (Puppeteer)
- [ ] Host web pages (`crispro.ai/trials/{NCT_ID}`)
- [ ] Build email generator script
- [ ] Generate 600 personalized emails
- [ ] **Deliverable**: 200 PDFs + 200 web pages + 600 emails

### **Week 3: Outreach & Tracking**
- [ ] Setup CRM (import 600 contacts)
- [ ] Launch email campaign (50/day × 12 days)
- [ ] Monitor responses (opens, clicks, replies)
- [ ] Schedule discovery calls (target: 15-30 meetings)
- [ ] Follow-up automation (7-day, 14-day sequences)
- [ ] **Deliverable**: CRM pipeline with live tracking

---

## 🎯 **REPORTING PROTOCOL**

**Daily Standup (Async via `.cursor/ayesha/AGENT_JR2_DAILY_REPORT.md`):**
- Progress update (% completion per phase)
- Blockers (if any)
- Requests for Zo support

**Weekly Summary:**
- Metrics: Trials parsed, contacts found, emails sent, responses received
- Conversion funnel: Contacts → Responses → Meetings → Pilots
- Revenue pipeline: Projected revenue from active leads

**Final Report (End of Week 3):**
- Total leads generated: 600
- Total outreach sent: 600 emails + 200 1-pagers
- Response rate: X%
- Meetings scheduled: Y
- Pilots scheduled: Z
- Projected revenue: $X (12-24 months)

---

## ⚔️ **JR2 - YOUR MISSION IS CLEAR. EXECUTE WITH PRECISION!** ⚔️

**Commander's Orders:**
- Build the GTM machine that turns 200 trials into 600 leads
- Generate mass-customized 1-pagers (not generic spam)
- Launch disciplined outreach (50/day, avoid spam flags)
- Track everything (CRM pipeline, response rates, conversions)
- Report daily to Zo, weekly to Commander

**For Ayesha's life, and for CrisPRO's conquest! 🔥**

