/**
 * Result types for analysis results
 */

import type { PaginatedResponse } from './common'

export type ResultType = 'qc_report' | 'alignment' | 'quantification' | 'de_analysis'

export type ChartType = 'line' | 'bar' | 'pie' | 'scatter' | 'histogram' | 'area'

export interface Result {
  id: string
  task_id: string
  result_type: ResultType
  result_path?: string
  metadata: Record<string, any>
  created_at: string
}

/**
 * @deprecated Use PaginatedResponse<Result> instead
 */
export type ResultListResponse = PaginatedResponse<Result>

// Chart data structures
export interface ChartData {
  type: ChartType
  x?: any[]
  y?: any[]
  categories?: string[]
  values?: number[]
  labels?: string[]
  bins?: number
  data?: number[]
}

// Visualization data
export interface QCMetrics {
  total_reads: number
  quality_score: number
  gc_content: number
  duplication_rate: number
}

export interface AlignmentStats {
  mapped_reads: number
  mapping_rate: number
  properly_paired: number
  average_coverage: number
}

export interface QuantificationMetrics {
  total_genes: number
  expressed_genes: number
  median_expression: number
}

export interface DEMetrics {
  total_genes: number
  up_regulated: number
  down_regulated: number
  significant: number
}

export interface GeneData {
  gene: string
  expression?: number
  log2_fold_change?: number
  p_value?: number
  significant?: boolean
}

export interface ResultVisualizationData {
  type: ResultType
  metrics: QCMetrics | AlignmentStats | QuantificationMetrics | DEMetrics | Record<string, any>
  charts: Record<string, ChartData>
  status?: 'pass' | 'warning' | 'fail'
  message?: string
  raw_metadata?: Record<string, any>
  top_genes?: GeneData[]
  significant_genes?: GeneData[]
}

export interface ResultSummary {
  task_id: string
  total_results: number
  result_types: ResultType[]
  results_by_type?: Record<string, {
    count: number
    latest: string
  }>
  summary?: string
}
