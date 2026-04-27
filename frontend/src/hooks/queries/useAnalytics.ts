/**
 * TanStack Query hooks for the Analytics + Knowledge Base domains.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import analyticsService from '@/services/analytics.service'

// ----- analytics ----------------------------------------------------------

export function useAnalyticsSummary(
  period: 'today' | 'week' | 'month' | 'year' | 'all' = 'week',
) {
  return useQuery({
    queryKey: [...queryKeys.analytics.dashboard, period],
    queryFn: () => analyticsService.getAnalyticsSummary(period),
    staleTime: 60_000,
    refetchInterval: 5 * 60_000,
  })
}

export function useTimeSeriesData(
  metric: string,
  startDate: string,
  endDate: string,
  interval: 'hour' | 'day' | 'week' | 'month' = 'day',
) {
  return useQuery({
    queryKey: queryKeys.analytics.timeSeries(metric, `${startDate}_${endDate}_${interval}`),
    queryFn: () => analyticsService.getTimeSeriesData(metric, startDate, endDate, interval),
    enabled: !!metric && !!startDate && !!endDate,
    staleTime: 60_000,
  })
}

export function useTrendAnalysis(metrics: string[]) {
  return useQuery({
    queryKey: ['analytics', 'trends', metrics],
    queryFn: () => analyticsService.getTrendAnalysis(metrics),
    enabled: metrics.length > 0,
    staleTime: 60_000,
  })
}

export function useProjectPerformance(projectId: string | undefined) {
  return useQuery({
    queryKey: queryKeys.analytics.projectPerformance(projectId ?? ''),
    queryFn: () => analyticsService.getProjectPerformance(projectId as string),
    enabled: !!projectId,
  })
}

export function useCompareProjects(projectIds: string[]) {
  return useQuery({
    queryKey: ['analytics', 'compare', projectIds],
    queryFn: () => analyticsService.compareProjects(projectIds),
    enabled: projectIds.length >= 2,
  })
}

export function useGenerateReport() {
  return useMutation({
    mutationFn: (config: any) => analyticsService.generateReport(config),
  })
}

// ----- knowledge base ------------------------------------------------------

export function useKnowledgeSearch(query: string, category?: any) {
  return useQuery({
    queryKey: ['knowledge', 'search', query, category],
    queryFn: () => analyticsService.searchKnowledge(query, category),
    enabled: query.length > 1,
    staleTime: 5 * 60_000,
  })
}

export function useArticlesByCategory(category: any) {
  return useQuery({
    queryKey: ['knowledge', 'category', category],
    queryFn: () => analyticsService.getArticlesByCategory(category),
    enabled: !!category,
    staleTime: 10 * 60_000,
  })
}

export function usePopularArticles(limit = 10) {
  return useQuery({
    queryKey: ['knowledge', 'popular', limit],
    queryFn: () => analyticsService.getPopularArticles(limit),
    staleTime: 30 * 60_000,
  })
}

export function useArticle(id: string | undefined) {
  return useQuery({
    queryKey: ['knowledge', 'article', id],
    queryFn: () => analyticsService.getArticle(id as string),
    enabled: !!id,
    staleTime: 10 * 60_000,
  })
}

export function useTutorials(category?: any) {
  return useQuery({
    queryKey: ['knowledge', 'tutorials', category],
    queryFn: () => analyticsService.getTutorials(category),
    staleTime: 10 * 60_000,
  })
}

export function useTutorialProgress() {
  return useQuery({
    queryKey: ['knowledge', 'tutorial-progress'],
    queryFn: () => analyticsService.getTutorialProgress(),
  })
}

export function useMarkArticleHelpful() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ articleId, helpful }: { articleId: string; helpful: boolean }) =>
      analyticsService.markArticleHelpful(articleId, helpful),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['knowledge'] })
    },
  })
}

export function useCompleteTutorialStep() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ tutorialId, stepId }: { tutorialId: string; stepId: string }) =>
      analyticsService.completeTutorialStep(tutorialId, stepId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['knowledge', 'tutorial-progress'] })
    },
  })
}
