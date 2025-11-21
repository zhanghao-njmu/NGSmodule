"""
Initialize built-in pipeline templates

This script creates default pipeline templates based on existing NGS analysis scripts.
"""
from app.core.database import SessionLocal
from app.models.pipeline_template import PipelineTemplate
import uuid


BUILTIN_TEMPLATES = [
    {
        "name": "preAlignmentQC",
        "display_name": "Pre-Alignment Quality Control",
        "description": "Perform quality control on raw sequencing data before alignment. Includes FastQC, Trimmomatic, and quality filtering.",
        "category": "Quality Control",
        "script_name": "preAlignmentQC",
        "script_path": "GeneralSteps/preAlignmentQC.sh",
        "default_params": {
            "threads": 8,
            "trim_quality": 20,
            "min_length": 36,
            "adapter_removal": True
        },
        "param_schema": {
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 64, "default": 8},
            "trim_quality": {"type": "integer", "label": "Trim Quality Threshold", "min": 0, "max": 40, "default": 20},
            "min_length": {"type": "integer", "label": "Minimum Read Length", "min": 20, "max": 200, "default": 36},
            "adapter_removal": {"type": "boolean", "label": "Remove Adapters", "default": True}
        },
        "estimated_time": "30min-2hours",
        "min_memory_gb": "4GB",
        "min_cpu_cores": "4",
        "sort_order": "1",
        "tags": ["quality-control", "preprocessing"]
    },
    {
        "name": "alignment",
        "display_name": "Sequence Alignment",
        "description": "Align sequencing reads to reference genome using STAR, HISAT2, or Bowtie2. Supports both RNA-seq and DNA-seq data.",
        "category": "Alignment",
        "script_name": "Alignment",
        "script_path": "GeneralSteps/Alignment.sh",
        "default_params": {
            "aligner": "STAR",
            "threads": 16,
            "genome": "hg38",
            "max_mismatches": 2,
            "output_format": "BAM"
        },
        "param_schema": {
            "aligner": {"type": "select", "label": "Aligner Tool", "options": ["STAR", "HISAT2", "Bowtie2"], "default": "STAR"},
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 64, "default": 16},
            "genome": {"type": "select", "label": "Reference Genome", "options": ["hg38", "hg19", "mm10", "mm39"], "default": "hg38"},
            "max_mismatches": {"type": "integer", "label": "Max Mismatches", "min": 0, "max": 10, "default": 2},
            "output_format": {"type": "select", "label": "Output Format", "options": ["BAM", "SAM"], "default": "BAM"}
        },
        "estimated_time": "2-6hours",
        "min_memory_gb": "32GB",
        "min_cpu_cores": "8",
        "sort_order": "2",
        "tags": ["alignment", "core"]
    },
    {
        "name": "postAlignmentQC",
        "display_name": "Post-Alignment Quality Control",
        "description": "Quality control after alignment. Includes RSeQC, Qualimap, and duplication analysis.",
        "category": "Quality Control",
        "script_name": "postAlignmentQC",
        "script_path": "GeneralSteps/postAlignmentQC.sh",
        "default_params": {
            "threads": 8,
            "check_duplication": True,
            "check_coverage": True
        },
        "param_schema": {
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 32, "default": 8},
            "check_duplication": {"type": "boolean", "label": "Check Duplication", "default": True},
            "check_coverage": {"type": "boolean", "label": "Check Coverage", "default": True}
        },
        "estimated_time": "1-3hours",
        "min_memory_gb": "16GB",
        "min_cpu_cores": "4",
        "sort_order": "3",
        "tags": ["quality-control", "postprocessing"]
    },
    {
        "name": "quantification",
        "display_name": "Gene Expression Quantification",
        "description": "Quantify gene expression levels from aligned reads. Supports featureCounts, HTSeq, and RSEM.",
        "category": "RNA-seq",
        "script_name": "Quantification",
        "script_path": "Analysis/Quantification/Quantification.sh",
        "default_params": {
            "method": "featureCounts",
            "threads": 8,
            "feature_type": "exon",
            "strand_specific": 0
        },
        "param_schema": {
            "method": {"type": "select", "label": "Quantification Method", "options": ["featureCounts", "HTSeq", "RSEM"], "default": "featureCounts"},
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 32, "default": 8},
            "feature_type": {"type": "select", "label": "Feature Type", "options": ["exon", "gene", "transcript"], "default": "exon"},
            "strand_specific": {"type": "select", "label": "Strand Specificity", "options": [{"value": 0, "label": "Unstranded"}, {"value": 1, "label": "Forward"}, {"value": 2, "label": "Reverse"}], "default": 0}
        },
        "estimated_time": "30min-2hours",
        "min_memory_gb": "8GB",
        "min_cpu_cores": "4",
        "sort_order": "4",
        "tags": ["quantification", "rna-seq"]
    },
    {
        "name": "differentialExpression",
        "display_name": "Differential Expression Analysis",
        "description": "Identify differentially expressed genes between conditions using DESeq2, edgeR, or limma.",
        "category": "RNA-seq",
        "script_name": "DifferentialExpression",
        "script_path": "Analysis/DifferentialExpression/DifferentialExpression.sh",
        "default_params": {
            "method": "DESeq2",
            "padj_cutoff": 0.05,
            "lfc_threshold": 1.0,
            "normalization": "TMM"
        },
        "param_schema": {
            "method": {"type": "select", "label": "DE Method", "options": ["DESeq2", "edgeR", "limma"], "default": "DESeq2"},
            "padj_cutoff": {"type": "number", "label": "Adjusted P-value Cutoff", "min": 0, "max": 1, "step": 0.01, "default": 0.05},
            "lfc_threshold": {"type": "number", "label": "Log2 Fold Change Threshold", "min": 0, "max": 10, "step": 0.1, "default": 1.0},
            "normalization": {"type": "select", "label": "Normalization Method", "options": ["TMM", "RLE", "upperquartile"], "default": "TMM"}
        },
        "estimated_time": "10min-1hour",
        "min_memory_gb": "8GB",
        "min_cpu_cores": "2",
        "sort_order": "5",
        "tags": ["differential-expression", "rna-seq", "analysis"]
    },
    {
        "name": "cnvAnalysis",
        "display_name": "Copy Number Variation Analysis",
        "description": "Detect copy number variations from sequencing data using HMMcopy, DNAcopy, or CNVkit.",
        "category": "DNA-seq",
        "script_name": "CNVanalysis",
        "script_path": "Analysis/CNV/CNVanalysis.sh",
        "default_params": {
            "method": "HMMcopy",
            "window_size": 1000,
            "threads": 8
        },
        "param_schema": {
            "method": {"type": "select", "label": "CNV Method", "options": ["HMMcopy", "DNAcopy", "CNVkit"], "default": "HMMcopy"},
            "window_size": {"type": "integer", "label": "Window Size (bp)", "min": 100, "max": 10000, "default": 1000},
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 32, "default": 8}
        },
        "estimated_time": "2-6hours",
        "min_memory_gb": "16GB",
        "min_cpu_cores": "4",
        "sort_order": "6",
        "tags": ["cnv", "dna-seq", "variant-calling"]
    },
    {
        "name": "gatkGermlineVariant",
        "display_name": "GATK Germline Variant Calling",
        "description": "Call germline short variants (SNPs and indels) using GATK HaplotypeCaller.",
        "category": "DNA-seq",
        "script_name": "GATK-germline-short-variant",
        "script_path": "Analysis/GATK/GATK-germline-short-variant.sh",
        "default_params": {
            "threads": 16,
            "min_base_quality": 20,
            "min_mapping_quality": 20,
            "call_confidence": 30.0
        },
        "param_schema": {
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 64, "default": 16},
            "min_base_quality": {"type": "integer", "label": "Min Base Quality", "min": 0, "max": 40, "default": 20},
            "min_mapping_quality": {"type": "integer", "label": "Min Mapping Quality", "min": 0, "max": 60, "default": 20},
            "call_confidence": {"type": "number", "label": "Call Confidence Threshold", "min": 0, "max": 100, "default": 30.0}
        },
        "estimated_time": "4-12hours",
        "min_memory_gb": "32GB",
        "min_cpu_cores": "8",
        "sort_order": "7",
        "tags": ["variant-calling", "dna-seq", "gatk", "germline"]
    },
    {
        "name": "gatkSomaticVariant",
        "display_name": "GATK Somatic Variant Calling",
        "description": "Call somatic mutations using GATK Mutect2 for tumor-normal paired samples.",
        "category": "DNA-seq",
        "script_name": "GATK-somatic-short-variant",
        "script_path": "Analysis/GATK/GATK-somatic-short-variant.sh",
        "default_params": {
            "threads": 16,
            "min_allele_frequency": 0.01,
            "tumor_lod": 5.3
        },
        "param_schema": {
            "threads": {"type": "integer", "label": "CPU Threads", "min": 1, "max": 64, "default": 16},
            "min_allele_frequency": {"type": "number", "label": "Min Allele Frequency", "min": 0, "max": 1, "step": 0.001, "default": 0.01},
            "tumor_lod": {"type": "number", "label": "Tumor LOD Threshold", "min": 0, "max": 20, "default": 5.3}
        },
        "estimated_time": "6-24hours",
        "min_memory_gb": "64GB",
        "min_cpu_cores": "16",
        "sort_order": "8",
        "tags": ["variant-calling", "dna-seq", "gatk", "somatic", "cancer"]
    }
]


def init_pipeline_templates():
    """
    Initialize built-in pipeline templates
    """
    db = SessionLocal()

    try:
        print("Initializing pipeline templates...")

        for template_data in BUILTIN_TEMPLATES:
            # Check if template already exists
            existing = db.query(PipelineTemplate).filter(
                PipelineTemplate.name == template_data["name"]
            ).first()

            if existing:
                print(f"  ℹ️  Template '{template_data['name']}' already exists, skipping...")
                continue

            # Create new template
            template = PipelineTemplate(
                id=uuid.uuid4(),
                **template_data,
                is_active=True,
                is_builtin=True
            )

            db.add(template)
            print(f"  ✅ Created template: {template_data['display_name']}")

        db.commit()
        print(f"\n✅ Pipeline templates initialized successfully! ({len(BUILTIN_TEMPLATES)} templates)")

    except Exception as e:
        print(f"❌ Error initializing templates: {e}")
        db.rollback()
    finally:
        db.close()


if __name__ == "__main__":
    init_pipeline_templates()
