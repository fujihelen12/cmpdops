# CompoundOps Assistant

## AI-Powered Experimental Feasibility and Assay Design Tool

CompoundOps Assistant is a cloud-hosted scientific workflow application that helps researchers evaluate compound availability, assess experimental feasibility, design dose-response assays, generate intermediate dilution strategies, and create assay-ready plate maps—all within a single guided workflow.

This project was developed as part of a Human-Centered Design course and demonstrates how domain expertise, AI-assisted software development, and cloud deployment can transform a complex scientific workflow into an intuitive web application.

---

## Live Application

Access the deployed application through Streamlit Cloud:

**[Insert Streamlit URL Here]**

Example:

```
https://your-app-name.streamlit.app
```

---

# Problem Statement

Scientists often rely on multiple disconnected tools when planning experiments:

- Inventory spreadsheets
- Concentration calculators
- SOP documents
- Dilution worksheets
- Plate mapping templates
- Manual calculations

This process can be time-consuming, inefficient, and prone to human error.

CompoundOps Assistant was designed to consolidate these tasks into a single workflow-driven application.

---

# Key Features

## Inventory Search

- Search inventory for compounds
- View stock concentration
- Review available volume
- Display solvent information
- Display storage conditions
- Review compound handling guidance

---

## Assay Designer

- Configure dose-response experiments
- Define top concentration
- Select dilution factor
- Choose curve points
- Calculate direct transfer volumes
- Generate solvent-normalized backfill volumes

---

## Intermediate Strategy Generator

- Create practical intermediate dilution recipes
- Reduce excessive direct transfers
- Support limited-transfer workflows
- Improve assay feasibility planning

---

## Dose Response Plate Mapping

- Generate plate-ready assay layouts
- Configure replicates
- Assign manual controls
- Visualize final plate organization

---

## SOP Guidance

- Recommend compound handling guidance
- Display storage considerations
- Provide solvent-specific recommendations
- Support experimental planning decisions

---

# Technology Stack

### Programming Language

- Python

### Front-End Framework

- Streamlit

### Data Processing

- Pandas

### Visualization

- Plotly

### Version Control

- GitHub

### Deployment

- Streamlit Cloud

### Development Accelerator

- Microsoft Copilot

---

# Human-Centered Design Process

This project was developed using Human-Centered Design principles.

## User Research

The workflow was based on real-world compound management and assay planning activities.

## Problem Definition

Scientists often need to switch between multiple disconnected tools when evaluating experimental feasibility.

## How Might We Statement

> How might we help scientists rapidly evaluate compound availability, design experiments, generate dilution strategies, and create plate layouts from a single workflow?

## Prototype Refinement

The application evolved through multiple iterations based on workflow analysis and usability improvements.

Key refinements included:

- Improved inventory selection logic
- SOP integration
- Aliquot selection workflow
- Intermediate concentration generation
- Enhanced plate map visualization
- Browser-based deployment

---

# Deployment Journey

The project began as a local prototype running on a single laptop.

To make the application accessible to teammates and instructors:

1. Project files were organized into a GitHub repository.
2. Dependencies were documented using `requirements.txt`.
3. The repository was connected to Streamlit Cloud.
4. The application was deployed as a browser-accessible web application.

The final product can now be accessed without:

- Installing Python
- Using Git
- Running code locally
- Installing additional software

Users simply open a URL and use the application.

---

# Repository Structure

```text
cmpdops/
│
├── app.py
├── inventory.csv
├── requirements.txt
├── README.md
│
├── DMSO_SOP.docx
├── Storage_SOP.docx
├── FreezeThaw_SOP.docx
└── Compound_Solvent_Handling_Feasibility_and_Plate_Storage_SOP.docx
```

---

# Example Workflow

1. Search compounds in inventory
2. Select desired compound/aliquot
3. Configure assay parameters
4. Generate dose-response design
5. Create intermediate dilution strategy
6. Generate plate map
7. Review SOP recommendations
8. Export assay-ready plan

---

# Future Enhancements

Potential future improvements include:

- Echo worklist generation
- Inventory consumption forecasting
- Solubility risk prediction
- Natural-language assay requests
- AI-assisted SOP guidance
- Automated feasibility recommendations
- Multi-plate study design

---

# Educational Disclaimer

This repository was developed as an academic project and demonstration application.

The publicly shared version should only contain:

- Demonstration data
- Non-confidential information
- Educational content

No proprietary or confidential scientific information should be uploaded to the public repository.

---

# Lessons Learned

This project demonstrated that:

- Domain expertise can drive software innovation
- Human-Centered Design improves usability
- AI-assisted development accelerates prototyping
- Scientific workflows can be transformed into intuitive web tools
- Cloud deployment significantly improves accessibility and collaboration

---

# Author

**Helen Martone**

Senior Associate Scientist – Compound Management

Human-Centered Design Course Project

2026

---

## Acknowledgments

Special thanks to the course instructors, teammates, and Microsoft Copilot for supporting the development, refinement, and deployment of this project.
