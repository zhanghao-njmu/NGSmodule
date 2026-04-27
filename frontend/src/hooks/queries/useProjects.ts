/**
 * TanStack Query hooks for the Projects domain.
 *
 * Project service is built on the generic CRUD factory, so we adapt the
 * factory method names (getAll/getById/create/update/delete) to the
 * domain-friendly names used here.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import projectService from '@/services/project.service'

export function useProjectList(params?: any) {
  return useQuery({
    queryKey: queryKeys.projects.list(params),
    queryFn: () => projectService.getAll(params),
  })
}

export function useProject(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.projects.detail(id ?? ''),
    queryFn: () => projectService.getById(id as string),
    enabled: !!id,
  })
}

export function useProjectDetails(id: string | undefined) {
  return useQuery({
    queryKey: [...queryKeys.projects.detail(id ?? ''), 'full'],
    queryFn: () => projectService.getProjectDetails(id as string),
    enabled: !!id,
  })
}

export function useGlobalProjectStats() {
  return useQuery({
    queryKey: queryKeys.stats.projects,
    queryFn: () => projectService.getStats(),
    staleTime: 60_000,
  })
}

export function useCreateProject() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => projectService.create(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.projects })
    },
  })
}

export function useUpdateProject() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, data }: { id: string; data: any }) => projectService.update(id, data),
    onSuccess: (_, { id }) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.detail(id) })
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.all })
    },
  })
}

export function useDeleteProject() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => projectService.delete(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.projects })
    },
  })
}

export function useArchiveProject() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => projectService.archiveProject(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.all })
    },
  })
}

export function useRestoreProject() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => projectService.restoreProject(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.projects.all })
    },
  })
}
