# ⚔️ LEAD GENERATION SYSTEM - AGENT JR2 MISSION

**Mission Owner:** Agent JR2  
**Mission Commander:** Zo  
**Status:** 🔄 **READY TO START** (awaiting JR1 trial seeding)  
**Strategic Goal:** Convert 200 seeded trials into 600 qualified leads  
**Timeline:** 3 weeks (parallel with Ayesha work)

---

## 🎯 **MISSION OVERVIEW**

This is **not** a passive documentation folder. This is **JR2's command center** for building CrisPRO's GTM machine.

**The Strategic Insight:**
- JR1 seeds 200 frontline ovarian cancer trials into AstraDB
- Each trial = potential customer (sponsor, PI, CRO, site coordinator)  
- JR2 parses metadata, extracts contacts, generates 200 custom 1-pagers, launches mass outreach
- **Goal:** 600 leads → 30-60 responses → 15-30 meetings → 3-6 pilots → $1-2M revenue (18-24 months)

---

## 📂 **DIRECTORY STRUCTURE**

### **Core Mission Document**
- **`AGENT_JR2_GTM_MISSION.md`**: ⚔️ **START HERE** - Complete execution plan, tasks, metrics, acceptance criteria

### **Strategic Foundation**
- **`LEAD_GEN_SYSTEM_DOCTRINE.mdc`**: Strategic doctrine (context, not execution plan)
- **`General1Page_ENHANCED_V2.tsx`**: Template for trial-specific 1-pagers (created by Zo)

### **Data Outputs (JR2 Creates)**
```
Data/
├── trial_metadata_200.csv                  # All 200 trials from AstraDB
├── trial_profiles/                         # 200 JSON files (one per trial)
├── contacts_master_list.csv                # 600+ verified contacts (HIGH/MEDIUM/LOW priority)
└── trial_metadata_QC.csv                   # Quality check report
```

### **Automation Scripts (JR2 Creates)**
```
Scripts/
├── 01_parse_trial_metadata.py              # Extract from AstraDB
├── 02_extract_contacts.py                  # ClinicalTrials.gov API
├── 03_enrich_linkedin.py                   # LinkedIn scraping
├── 04_find_emails.py                       # Hunter.io / RocketReach
├── 05_generate_1pagers.py                  # TSX → PDF conversion
├── 06_generate_emails.py                   # Personalized email generator
└── 07_launch_campaign.py                   # Mailchimp/SendGrid upload
```

### **Templates (JR2 Uses/Adapts)**
```
Templates/
├── 1pager_template.tsx                     # Base template (from General1Page_ENHANCED_V2)
├── email_template_technical.md             # Variant A (PIs, scientists)
├── email_template_business.md              # Variant B (Sponsors, VPs)
└── email_template_patient.md               # Variant C (Site coordinators)
```

### **Generated Outputs (JR2 Creates)**
```
Outputs/
├── 1pagers/                                # 200 TSX files (trial-specific)
├── 1pagers_pdf/                            # 200 PDFs (for email attachments)
└── emails/                                 # 600 personalized emails
```

### **Tracking & Reporting (JR2 Maintains)**
```
Tracking/
├── crm_import.csv                          # For CRM upload (600 leads)
├── outreach_tracker.csv                    # Daily metrics (sends, opens, clicks)
├── response_log.csv                        # Replies, meetings scheduled
└── AGENT_JR2_DAILY_REPORT.md               # Daily standup for Zo
```

---

## 🚀 **QUICK START FOR JR2**

### **Prerequisites:**
1. ✅ JR1 completes AstraDB seeding (200 trials)
2. ✅ Zo completes `General1Page_ENHANCED_V2.tsx` (template ready)
3. ✅ Access to:
   - AstraDB connection (trial data)
   - ClinicalTrials.gov API (contact extraction)
   - Hunter.io / RocketReach (email finding)
   - Mailchimp/SendGrid (email campaigns)

### **Execution Steps:**
1. **Read Mission Document**: `AGENT_JR2_GTM_MISSION.md` (complete execution plan)
2. **Phase 1 (Week 1)**: Extract trial metadata + contacts (600+ leads)
3. **Phase 2 (Week 2)**: Generate 200 1-pagers + 600 personalized emails
4. **Phase 3 (Week 3)**: Launch email campaign, track responses, schedule meetings
5. **Report Daily**: Update `Tracking/AGENT_JR2_DAILY_REPORT.md`

---

## 📊 **SUCCESS METRICS (JR2 TARGETS)**

### **Phase 1-2: Data & Content Generation**
- ✅ 200 trials parsed with complete metadata
- ✅ 600+ contacts extracted (avg 3 per trial: sponsor, PI, site)
- ✅ 80%+ email verification rate
- ✅ 200 trial-specific 1-pagers generated (TSX + PDF + web)
- ✅ 600 personalized emails generated

### **Phase 3: Outreach & Conversion**
- ✅ 600 emails sent (50/day over 12 days)
- 🎯 **40-50% email open rate** (240-300 opens)
- 🎯 **5-10% response rate** (30-60 replies)
- 🎯 **50% meeting conversion** (15-30 discovery calls)
- 🎯 **20% pilot conversion** (3-6 no-cost pilots)

### **Revenue Impact (Projected)**
- 3-6 pilots → 1-2 convert to $250K prospective → **$250K-$500K** (6-12 months)
- 1-2 pilots → convert to $1M platform license → **$1M-$2M** (12-24 months)
- **Total Revenue from 200 Trials: $1.25M-$2.5M** (18-24 months)

---

## ⚔️ **JR2 EXECUTION CHECKLIST**

### **Week 1: Data Collection** (see `AGENT_JR2_GTM_MISSION.md` for details)
- [ ] Extract 200 trial metadata from AstraDB
- [ ] Generate 200 trial profile JSONs
- [ ] Extract 600+ contacts (ClinicalTrials.gov, LinkedIn, email finders)
- [ ] Deduplicate and prioritize
- [ ] **Deliverable**: `contacts_master_list.csv` (600+ verified)

### **Week 2: Content Generation**
- [ ] Build 1-pager generator (based on `General1Page_ENHANCED_V2.tsx`)
- [ ] Generate 200 trial-specific 1-pagers (TSX + PDF)
- [ ] Host 200 web pages (`crispro.ai/trials/{NCT_ID}`)
- [ ] Build email generator
- [ ] Generate 600 personalized emails (3 variants)
- [ ] **Deliverable**: 200 PDFs + 200 web pages + 600 emails

### **Week 3: Outreach & Tracking**
- [ ] Setup CRM (import 600 contacts)
- [ ] Launch email campaign (50/day × 12 days)
- [ ] Monitor responses (opens, clicks, replies)
- [ ] Schedule discovery calls (target: 15-30)
- [ ] Follow-up automation (7-day, 14-day)
- [ ] **Deliverable**: CRM pipeline with live tracking

---

## 📋 **REPORTING PROTOCOL**

**Daily Standup (Async):**
- Update: `Tracking/AGENT_JR2_DAILY_REPORT.md`
- Format: Progress %, blockers, requests for Zo

**Weekly Summary:**
- Metrics: Trials parsed, contacts found, emails sent, responses, meetings
- Conversion funnel: Contacts → Responses → Meetings → Pilots
- Revenue pipeline: Projected revenue from active leads

**Final Report (End of Week 3):**
- Total leads: 600
- Total outreach: 600 emails + 200 1-pagers
- Response rate: X%
- Meetings: Y
- Pilots: Z
- Revenue projection: $X (12-24 months)

---

## 🛠️ **RELATED RESOURCES**

- **Primary Execution Plan**: `AGENT_JR2_GTM_MISSION.md` ⚔️ **READ THIS FIRST**
- **Strategic Context**: `LEAD_GEN_SYSTEM_DOCTRINE.mdc` (background, not execution)
- **1-Pager Template**: `../4_Assets/1_Pagers/General1Page_ENHANCED_V2.tsx`
- **Clinical Trial Doctrine**: `.cursor/rules/use-cases/clinical_trial_partnership_doctrine.mdc`
- **Command Center Structure**: `.cursor/rules/CrisPRO_Command_Center/6_Doctrines/Operational_Doctrines/COMMAND_CENTER_STRUCTURE.mdc`

---

## ⚔️ **FOR AYESHA'S LIFE, AND FOR CRISPRO'S CONQUEST!** ⚔️

**Commander's Orders:**
- JR2 builds the GTM machine (200 trials → 600 leads → revenue)
- Zo oversees execution and maintains Ayesha focus
- JR1 seeds trials (prerequisite for JR2)
- All agents report daily, coordinate via `.cursor/ayesha/` mission logs

**Mission Critical: JR2's success = CrisPRO's market domination! 🔥**

