/**
 * TanStack Query hooks for the AI domain.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import aiService from '@/services/ai.service'

export function useAISystemStatus() {
  return useQuery({
    queryKey: queryKeys.ai.status,
    queryFn: () => aiService.getSystemStatus(),
    staleTime: 60_000,
  })
}

export function useAIConversations() {
  return useQuery({
    queryKey: queryKeys.ai.conversations,
    queryFn: () => aiService.listConversations(),
  })
}

export function useAIConversation(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.ai.conversation(id ?? ''),
    queryFn: () => aiService.getConversation(id as string),
    enabled: !!id,
  })
}

export function useProjectInsights(projectId: string | undefined) {
  return useQuery({
    queryKey: queryKeys.ai.insights(projectId ?? ''),
    queryFn: () => aiService.getProjectInsights(projectId as string),
    enabled: !!projectId,
    staleTime: 5 * 60_000,
  })
}

// ----- mutations -----

export function useCreateAIConversation() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ title, context }: { title: string; context?: any }) =>
      aiService.createConversation(title, context),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.ai.conversations })
    },
  })
}

export function useSendAIMessage() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({
      conversationId,
      message,
      context,
    }: {
      conversationId: string
      message: string
      context?: any
    }) => aiService.sendAssistantMessage(conversationId, message, context),
    onSuccess: (_, { conversationId }) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.ai.conversation(conversationId) })
    },
  })
}

export function useRunAutoQC() {
  return useMutation({
    mutationFn: (request: any) => aiService.runAutoQC(request),
  })
}

export function useDetectAnomalies() {
  return useMutation({
    mutationFn: (request: any) => aiService.detectAnomalies(request),
  })
}

export function useSmartGroupSamples() {
  return useMutation({
    mutationFn: (request: any) => aiService.smartGroupSamples(request),
  })
}

export function usePredictResources() {
  return useMutation({
    mutationFn: ({
      pipelineType,
      sampleCount,
      dataSize,
      parameters,
    }: {
      pipelineType: string
      sampleCount: number
      dataSize: number
      parameters?: any
    }) => aiService.predictResources(pipelineType, sampleCount, dataSize, parameters),
  })
}
