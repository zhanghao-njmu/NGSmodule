/**
 * TanStack Query hooks for the Pipelines domain (templates + execution).
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import pipelineService from '@/services/pipeline.service'

export function usePipelineTemplates(params?: any) {
  return useQuery({
    queryKey: [...queryKeys.pipelines.templates, params],
    queryFn: () => pipelineService.getAll(params),
    // Templates are mostly static — cache aggressively
    staleTime: 5 * 60_000,
  })
}

export function usePipelineTemplate(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.pipelines.template(id ?? ''),
    queryFn: () => pipelineService.getById(id as string),
    enabled: !!id,
    staleTime: 5 * 60_000,
  })
}

export function usePipelineCategories() {
  return useQuery({
    queryKey: [...queryKeys.pipelines.all, 'categories'],
    queryFn: () => pipelineService.getCategories(),
    staleTime: 10 * 60_000,
  })
}

// ----- mutations -----

export function useExecutePipeline() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => pipelineService.executePipeline(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.tasks })
    },
  })
}

export function useBatchExecutePipeline() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => pipelineService.batchExecutePipeline(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.tasks })
    },
  })
}

export function usePipelineParameterRecommendations(templateId: string | undefined, params?: any) {
  return useQuery({
    queryKey: [...queryKeys.pipelines.all, 'recommendations', templateId, params],
    queryFn: () => pipelineService.getParameterRecommendations(templateId as string, params),
    enabled: !!templateId,
  })
}
